"""Thickness / layer-number maps from an optical micrograph and a simulated TCD reference.

Steps
  1. Read the sRGB micrograph, optionally downsample (box filter) and smooth (Gaussian).
  2. Linearise (sRGB EOTF) and calibrate the camera against the simulation: the bare
     substrate region (and optional extra regions of known structure) must reproduce
     the simulated colours. One anchor -> per-channel gains (white balance + exposure);
     three or more anchors -> a full 3x3 colour-correction matrix.
  3. Convert every pixel to CAM02-UCS.
  4. TCD map = Delta E (CAM02-UCS) of each pixel from the simulated substrate, in
     absolute units (no per-image max normalisation).
  5. Assign each pixel the reference entry nearest in the full J'a'b' space, which
     separates structures with equal Delta E but different hue. Pixels whose nearest
     reference is further than ``max_residual`` (dust, edges, saturated pixels, colours the
     model cannot produce) are left unassigned. The best match whose label differs by more
     than the reference's ambiguity gap is also reported, since interference colours repeat
     with thickness: where it fits almost as well, the colour alone cannot decide.
  6. For references whose label is a whole number of layers (``label.classes``), the label is
     rounded half up into classes and cleaned with a majority filter.
  7. Every assigned pixel gets a reliability score (fit x uniqueness x consistency, see
     reliability.py). Optional checks, regions of known thickness (AFM), are compared with the
     map and set the colour error that the score assumes.

Nothing here depends on the material system: everything comes from the reference.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

from . import colorimetry as C
from . import reliability as REL
from .reference import Reference

Roi = tuple[int, int, int, int]  # x, y, width, height in pixels of the original image


@dataclass
class Anchor:
    roi: Roi
    ref_index: int
    name: str


@dataclass
class TCDResult:
    image: np.ndarray
    calibrated: np.ndarray
    Jab: np.ndarray
    tcd: np.ndarray
    index: np.ndarray
    value: np.ndarray
    residual: np.ndarray
    alt_value: np.ndarray | None
    alt_residual: np.ndarray | None
    saturated: np.ndarray
    reference: Reference
    anchors: list[Anchor]
    correction: np.ndarray
    scale: float
    settings: dict = field(default_factory=dict)
    reliability: np.ndarray | None = None   # 0-1 score, NaN where unassigned
    fit: np.ndarray | None = None
    uniqueness: np.ndarray | None = None
    consistency: np.ndarray | None = None
    checks: list[dict] = field(default_factory=list)
    sigma: dict = field(default_factory=dict)  # noise, model, total colour error (Delta E per axis)

    @property
    def assigned(self) -> np.ndarray:
        return self.index >= 0


def load_image(path: str | Path, max_side: int | None = 1600) -> tuple[np.ndarray, float]:
    """sRGB image as float in [0, 1] and the downsampling factor applied (<= 1)."""
    im = Image.open(path)
    if im.mode in ("I;16", "I;16B", "I"):
        arr = np.asarray(im, dtype=float)
        arr = np.repeat((arr / arr.max())[..., None], 3, axis=2)
        im = Image.fromarray((arr * 255).astype(np.uint8))
    im = im.convert("RGB")
    scale = 1.0
    if max_side and max(im.size) > max_side:
        scale = max_side / max(im.size)
        im = im.resize((round(im.width * scale), round(im.height * scale)), Image.BOX)
    return np.asarray(im, dtype=float) / 255.0, scale


def _roi_pixels(img: np.ndarray, roi: Roi, scale: float) -> np.ndarray:
    x, y, w, h = (int(round(v * scale)) for v in roi)
    patch = img[max(y, 0):y + max(h, 1), max(x, 0):x + max(w, 1)]
    if patch.size == 0:
        raise ValueError(f"ROI {roi} lies outside the image")
    return patch.reshape(-1, 3)


def fit_correction(measured: np.ndarray, target: np.ndarray, mode: str = "auto",
                   ridge: float = 1e-2) -> np.ndarray:
    """3x3 matrix M with target ~ measured @ M.T (linear RGB).

    'diagonal' : per-channel gains (least squares over the anchors).
    'matrix'   : full matrix, ridge-regularised towards the diagonal solution, so two
                 anchors (e.g. substrate + a known monolayer) already give a well-posed fit.
    'auto'     : diagonal for one anchor, matrix otherwise.
    """
    measured, target = np.atleast_2d(measured), np.atleast_2d(target)
    D = np.diag(np.sum(measured * target, axis=0) / np.sum(measured**2, axis=0))
    if mode == "diagonal" or (mode == "auto" and len(measured) == 1):
        return D
    G = measured.T @ measured
    lam = ridge * np.trace(G) / 3
    return np.linalg.solve(G + lam * np.eye(3), measured.T @ target + lam * D.T).T


def mode_filter(labels: np.ndarray, size: int, n_classes: int) -> np.ndarray:
    """Majority filter for an integer label map (-1 = unassigned is kept as its own class)."""
    shifted = labels + 1
    counts = np.stack([ndimage.uniform_filter((shifted == c).astype(float), size)
                       for c in range(n_classes + 1)])
    return counts.argmax(axis=0) - 1


def analyse(image_path: str | Path, reference: Reference, substrate_roi: Roi,
            anchors: list[tuple[Roi, int, str]] | None = None, max_side: int | None = 1600,
            smooth_sigma: float = 1.2, max_residual: float = 15.0, correction_mode: str = "auto",
            ambiguity_gap: float | None = None, saturation_level: float = 0.985,
            gamma: float | str = 1.0, checks: list[tuple[Roi, float, str]] | None = None,
            colour_error: float | None = None, tolerance: float | None = None) -> TCDResult:
    """Map ``image_path`` with ``reference``; the substrate region calibrates the camera.

    ``gamma`` decodes camera values to linear intensity: 'srgb' for the sRGB curve, or a
    number g for v**g. An AFM check of MoO3 flakes confirmed g = 1 for the camera used so far.
    ``ambiguity_gap`` (label units) defaults to the reference's ``label['gap']``.
    ``checks``: (region, known label value, name) of regions measured another way (AFM). They
    are not used for calibration; they are compared with the map.
    ``colour_error``: model colour error (Delta E per axis) for the reliability score; default
    is the estimate from the checks, or 3.0 without checks. ``tolerance``: how close (label
    units) counts as right; default half the gap, or the same class.
    """
    gap = float(reference.label["gap"] if ambiguity_gap is None else ambiguity_gap)
    tol = float(REL.default_tolerance(reference) if tolerance is None else tolerance)
    img, scale = load_image(image_path, max_side)
    saturated = np.any(img >= saturation_level, axis=-1)
    work = ndimage.gaussian_filter(img, sigma=(smooth_sigma, smooth_sigma, 0)) if smooth_sigma > 0 else img
    lin = C.srgb_to_linear(work) if gamma == "srgb" else np.power(work, float(gamma))

    all_anchors = [Anchor(tuple(substrate_roi), 0, "substrate")]
    all_anchors += [Anchor(tuple(r), i, n) for r, i, n in (anchors or [])]
    measured = np.array([np.median(_roi_pixels(lin, a.roi, scale), axis=0) for a in all_anchors])
    target = C.XYZ_to_linear(reference.XYZ[[a.ref_index for a in all_anchors]])
    M = fit_correction(measured, target, correction_mode)

    lin_cal = np.clip(lin @ M.T, 0, None)
    Jab = C.XYZ_to_cam02ucs(C.linear_to_XYZ(lin_cal))
    Jab = np.nan_to_num(Jab, nan=0.0)
    tcd = C.delta_e(Jab, reference.Jab[0])

    flat = Jab.reshape(-1, 3)
    tree = cKDTree(reference.Jab)
    residual, index = tree.query(flat)
    value = reference.value[index].astype(float)
    bad = (residual > max_residual) | saturated.ravel()
    index = np.where(bad, -1, index)
    value[bad] = np.nan
    value_map, residual_map = value.reshape(tcd.shape), residual.reshape(tcd.shape)

    # colour error for the reliability score: camera noise (substrate) + model error (checks or default)
    noise = REL.noise_sigma(Jab, _roi_pixels(Jab, substrate_roi, scale))
    check_rows = []
    for roi, known, name in checks or []:
        x, y, w, h = (int(round(v * scale)) for v in roi)
        sl = (slice(max(y, 0), y + max(h, 1)), slice(max(x, 0), x + max(w, 1)))
        check_rows.append((name, float(known), sl))
    nan_map = np.full(tcd.shape, np.nan)
    prelim = [REL.check_statistics(n, k, value_map[sl], residual_map[sl], nan_map[sl], reference, tol)
              for n, k, sl in check_rows]
    estimated = REL.model_error_from_checks(prelim, noise)
    model = float(colour_error if colour_error is not None else estimated if estimated is not None
                  else REL.DEFAULT_MODEL_ERROR)
    sigma = float(np.hypot(noise, model))

    alt_v, alt_d, uniq = REL.match_statistics(flat, index, ~bad, reference, gap, tol, sigma)
    alt_value, alt_residual = alt_v.reshape(tcd.shape), alt_d.reshape(tcd.shape)
    fit = np.where(np.isfinite(value_map), REL.chi3_sf(residual_map / sigma), np.nan)
    uniq = uniq.reshape(tcd.shape)
    cons = REL.consistency(value_map, reference, tol)
    score = fit * uniq * cons

    checks_out = []
    for n, k, sl in check_rows:
        row = REL.check_statistics(n, k, value_map[sl], residual_map[sl], score[sl], reference, tol)
        row.pop("right_residuals")
        checks_out.append(dict(row, roi=next(r for r, kk, nn in checks if nn == n)))

    calibrated_display = np.clip(C.colour.models.eotf_inverse_sRGB(np.clip(lin_cal, 0, 1)), 0, 1)
    settings = dict(max_side=max_side, smooth_sigma=smooth_sigma, max_residual=max_residual,
                    correction_mode=correction_mode, ambiguity_gap=gap, gamma=gamma, tolerance=tol,
                    image=str(image_path))
    sig = dict(noise=noise, model=model, total=sigma,
               model_source="given" if colour_error is not None else "checks" if estimated is not None else "default")
    return TCDResult(img, calibrated_display, Jab, tcd, index.reshape(tcd.shape), value_map, residual_map,
                     alt_value, alt_residual, saturated, reference, all_anchors, M, scale, settings,
                     score, fit, uniq, cons, checks_out, sig)


def layer_classes(result: TCDResult, clean: int = 5) -> np.ndarray:
    """Whole-class map for ``label.classes`` references (-1 = unassigned), majority-filtered.

    The label is rounded half up (1.5 -> 2) so both languages agree; classes start at the
    smallest class of the reference (normally 0 = substrate).
    """
    ref = result.reference
    lo = ref.class_range().start
    cls = np.where(result.assigned, np.floor(np.nan_to_num(result.value, nan=0) + 0.5).astype(int) - lo, -1)
    n = len(ref.class_range())
    cls = mode_filter(cls, clean, n) if clean and clean > 1 else cls
    return np.where(cls >= 0, cls + lo, -1)


def nice_step(span: float, n: int = 12) -> float:
    """A round bin width (1, 2, 2.5 or 5 x 10^k) giving about ``n`` bins over ``span``."""
    raw = max(span, 1e-12) / n
    p = 10 ** np.floor(np.log10(raw))
    return float(next(m * p for m in (1, 2, 2.5, 5, 10) if m * p >= raw - 1e-12))


def summary_rows(result: TCDResult) -> list[dict]:
    """Pixel fractions per class or per label bin, plus how often the colour is ambiguous.

    'covered' pixels are assigned pixels whose label differs from the calibration row's by
    more than half the ambiguity gap (the flake / bead area). Of those, the ambiguity row
    counts the ones whose best alternative fits within 2 Delta E of the chosen match.
    """
    ref, v = result.reference, result.value
    total = result.index.size
    gap = result.settings["ambiguity_gap"]

    def stats(item, m, of=total):
        return dict(item=item, pixels=int(m.sum()), fraction=float(m.sum() / max(of, 1)),
                    median_residual=float(np.median(result.residual[m])) if m.any() else np.nan)

    rows = [dict(item="unassigned", pixels=int((~result.assigned).sum()),
                 fraction=float((~result.assigned).mean()), median_residual=np.nan)]
    covered = result.assigned & (np.abs(np.nan_to_num(v, nan=ref.value[0]) - ref.value[0]) > gap / 2)
    close = covered & (result.alt_residual - result.residual < 2.0)
    unit = ref.label.get("unit", "")
    rows.append(dict(stats(f"covered px with an alternative > {gap:g} {unit} away within 2 dE".replace("  ", " "),
                           close, int(covered.sum())), median_residual=np.nan))
    if result.reliability is not None:
        reliable = covered & (np.nan_to_num(result.reliability) >= 0.5)
        rows.append(dict(stats("covered px with reliability >= 0.5", reliable, int(covered.sum())),
                         median_residual=np.nan))
    if ref.classes:
        cls = layer_classes(result, clean=0)
        for c in ref.class_range():
            rows.append(stats(ref.class_name(c), cls == c))
    else:
        w = nice_step(float(np.nanmax(ref.value) - np.nanmin(ref.value)))
        edges = np.arange(np.floor(np.nanmin(ref.value) / w) * w, np.nanmax(ref.value) + w * 0.999, w)
        for j, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            m = (v >= lo) & ((v < hi) if j < len(edges) - 2 else (v <= hi))
            rows.append(stats(f"{lo:g}-{hi:g} {unit}".strip(), m))
    return rows
