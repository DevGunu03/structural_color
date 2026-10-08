"""How much to trust each pixel of a thickness / layer map, and checks against known regions.

Reliability score (0-1) = fit x uniqueness x consistency, per assigned pixel:

  fit          Could the matched structure have produced this colour? The probability of a
               colour error at least as large as the residual, if errors are Gaussian with
               standard deviation sigma per J'a'b' axis (the chi distribution with 3 degrees of
               freedom): fit = erfc(x / sqrt 2) + sqrt(2 / pi) x exp(-x^2 / 2), x = residual / sigma.
  uniqueness   Does only the reported value explain the colour? Every candidate k gets the
               weight w_k = p_k exp(-d_k^2 / 2 sigma^2) (d_k = colour distance, p_k = a prior that
               is uniform in the label). uniqueness = the share of the weight that lies within
               the tolerance of the reported value (or in the same class).
  consistency  Do the neighbours agree? The share of assigned pixels in the 5 x 5 window whose
               value is within the tolerance (or in the same class).

sigma combines the camera noise, measured in the substrate region, and the colour error of the
model itself (wrong optical constants, oxide thickness, NA, ...): sigma^2 = noise^2 + model^2.
The model error defaults to 3 Delta E: on MoO3 flakes measured by AFM, correctly mapped pixels
were a median 4.4 Delta E from their simulated colour (= 1.54 sigma for 3-D Gaussian errors).
Checks (regions of known thickness) measure it for your own images instead.

The score is computed from the model. It cannot see an error of the model itself: if the
optical constants or the oxide thickness are wrong, a pixel can match a wrong thickness that
fits well and has no look-alike, and still score high. Checks are the only protection.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage, special

DEFAULT_MODEL_ERROR = 3.0          # Delta E (CAM02-UCS), see module docstring
CHI3_MEDIAN = 1.5381722            # median of the chi distribution with 3 degrees of freedom


def chi3_sf(x: np.ndarray) -> np.ndarray:
    """P(|e| >= x sigma) for a 3-D Gaussian error with sigma per axis."""
    x = np.asarray(x, dtype=float)
    return special.erfc(x / np.sqrt(2)) + np.sqrt(2 / np.pi) * x * np.exp(-x ** 2 / 2)


def label_prior(values: np.ndarray) -> np.ndarray:
    """Candidate weights that are uniform in the label, whatever the spacing of the sweep.

    Each row gets half the distance to each neighbour (in label order), each half capped at
    twice the typical step, so an isolated row (e.g. a bare substrate before a gap) does not
    get a huge share.
    """
    v = np.asarray(values, dtype=float)
    order = np.argsort(v, kind="stable")
    s = v[order]
    gaps = np.diff(s)
    typical = np.median(gaps[gaps > 0]) if np.any(gaps > 0) else 1.0
    left = np.minimum(np.r_[typical, gaps], 2 * typical)
    right = np.minimum(np.r_[gaps, typical], 2 * typical)
    w = np.empty_like(s)
    w[order] = np.maximum((left + right) / 2, 1e-12 * typical)
    return w / w.sum()


def same(ref, a: np.ndarray, b: np.ndarray, tol: float) -> np.ndarray:
    """True where labels a and b count as the same structure (same class, or within tol)."""
    if ref.classes:
        return np.floor(a + 0.5) == np.floor(b + 0.5)
    return np.abs(a - b) <= tol


def default_tolerance(ref) -> float:
    """Half the ambiguity gap (e.g. +/-20 nm for MoO3); for classes the class itself."""
    return 0.5 if ref.classes else ref.label["gap"] / 2


def match_statistics(flat: np.ndarray, index: np.ndarray, valid: np.ndarray, ref, gap: float, tol: float,
                     sigma: float, chunk: int = 20000) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Best alternative more than ``gap`` away (value, distance) and uniqueness, for each valid pixel.

    Exhaustive over all candidates: with fine thickness steps the nearest neighbours in colour
    all sit next to the best match, so a k-nearest search would miss the look-alikes.
    """
    n = flat.shape[0]
    alt_v, alt_d, uniq = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    R, vals = ref.Jab, ref.value
    rr = np.sum(R ** 2, axis=1)
    prior = label_prior(vals)
    idx = np.flatnonzero(valid)
    for s in range(0, len(idx), chunk):
        sel = idx[s:s + chunk]
        x = flat[sel]
        d2 = np.maximum(np.sum(x ** 2, axis=1)[:, None] + rr[None, :] - 2 * x @ R.T, 0)
        best = vals[index[sel]][:, None]
        w = prior[None, :] * np.exp(-(d2 - d2.min(axis=1, keepdims=True)) / (2 * sigma ** 2))
        uniq[sel] = np.sum(w * same(ref, vals[None, :], best, tol), axis=1) / np.sum(w, axis=1)
        d2[np.abs(vals[None, :] - best) <= gap] = np.inf
        j = np.argmin(d2, axis=1)
        dmin = d2[np.arange(len(sel)), j]
        ok = np.isfinite(dmin)
        alt_v[sel[ok]] = vals[j[ok]]
        alt_d[sel[ok]] = np.sqrt(dmin[ok])
    return alt_v, alt_d, uniq


def consistency(value: np.ndarray, ref, tol: float, size: int = 5) -> np.ndarray:
    """Share of assigned pixels in the size x size window that agree with the centre pixel."""
    H, W = value.shape
    r = size // 2
    padded = np.pad(value, r, constant_values=np.nan)
    agree = np.zeros((H, W))
    count = np.zeros((H, W))
    for dy in range(size):
        for dx in range(size):
            nb = padded[dy:dy + H, dx:dx + W]
            ok = np.isfinite(nb)
            count += ok
            with np.errstate(invalid="ignore"):
                agree += ok & same(ref, nb, value, tol)
    return np.where(np.isfinite(value), agree / np.maximum(count, 1), np.nan)


def noise_sigma(Jab: np.ndarray, roi_pixels: np.ndarray) -> float:
    """Per-axis colour noise (Delta E) of the substrate region: RMS scatter / sqrt 3."""
    d = roi_pixels - np.median(roi_pixels, axis=0)
    return float(np.sqrt(np.mean(np.sum(d ** 2, axis=1)) / 3))


def check_statistics(name: str, known: float, value: np.ndarray, residual: np.ndarray, reliability: np.ndarray,
                     ref, tol: float) -> dict:
    """How the map compares with a region of known structure."""
    assigned = np.isfinite(value)
    v, res, rel = value[assigned], residual[assigned], reliability[assigned]
    right = same(ref, v, np.full_like(v, known), tol)
    return dict(check=name, known=known, pixels=int(value.size), assigned=int(assigned.sum()),
                median_value=float(np.median(v)) if v.size else np.nan,
                median_abs_error=float(np.median(np.abs(v - known))) if v.size else np.nan,
                within_tolerance=float(right.mean()) if v.size else np.nan,
                median_reliability=float(np.median(rel)) if v.size else np.nan,
                median_residual=float(np.median(res)) if v.size else np.nan,
                right_residuals=res[right])


def model_error_from_checks(checks: list[dict], noise: float, min_pixels: int = 30) -> float | None:
    """Colour error of the model from check pixels that were mapped to the right structure.

    For 3-D Gaussian errors the median distance is 1.54 sigma, so sigma_total = median / 1.54
    and model = sqrt(sigma_total^2 - noise^2). None if too few check pixels were right.
    """
    r = np.concatenate([c["right_residuals"] for c in checks]) if checks else np.array([])
    if r.size < min_pixels:
        return None
    total = np.median(r) / CHI3_MEDIAN
    return float(np.sqrt(max(total ** 2 - noise ** 2, 0.25)))
