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
     model cannot produce) are left unassigned. For MoO3 the best match whose thickness
     differs by more than ``ambiguity_gap_nm`` is also reported, since interference
     colours repeat with thickness.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

from . import colorimetry as C
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


def best_alternative(flat: np.ndarray, index: np.ndarray, valid: np.ndarray, reference: Reference,
                     gap: float, chunk: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    """For each valid pixel, the closest reference row whose value differs from the best match
    by more than ``gap`` (exhaustive over all rows: with fine thickness steps the nearest
    neighbours all sit next to the best match, so a k-nearest search misses the look-alikes)."""
    alt_v = np.full(flat.shape[0], np.nan)
    alt_d = np.full(flat.shape[0], np.nan)
    ref, vals = reference.Jab, reference.value
    rr = np.sum(ref**2, axis=1)
    idx = np.flatnonzero(valid)
    for s in range(0, len(idx), chunk):
        sel = idx[s:s + chunk]
        x = flat[sel]
        d2 = np.sum(x**2, axis=1)[:, None] + rr[None, :] - 2 * x @ ref.T
        near = np.abs(vals[None, :] - vals[index[sel]][:, None]) <= gap
        d2[near] = np.inf
        j = np.argmin(d2, axis=1)
        dmin = d2[np.arange(len(sel)), j]
        ok = np.isfinite(dmin)
        alt_v[sel[ok]] = vals[j[ok]]
        alt_d[sel[ok]] = np.sqrt(np.maximum(dmin[ok], 0))
    return alt_v, alt_d


def mode_filter(labels: np.ndarray, size: int, n_classes: int) -> np.ndarray:
    """Majority filter for an integer label map (-1 = unassigned is kept as its own class)."""
    shifted = labels + 1
    counts = np.stack([ndimage.uniform_filter((shifted == c).astype(float), size)
                       for c in range(n_classes + 1)])
    return counts.argmax(axis=0) - 1


def analyse(image_path: str | Path, reference: Reference, substrate_roi: Roi,
            anchors: list[tuple[Roi, int, str]] | None = None, max_side: int | None = 1600,
            smooth_sigma: float = 1.2, max_residual: float = 15.0, correction_mode: str = "auto",
            ambiguity_gap_nm: float = 40.0, saturation_level: float = 0.985,
            gamma: float | str = 1.0) -> TCDResult:
    """``gamma`` decodes camera values to linear intensity: 'srgb' for the sRGB curve, or a
    number g for v**g. Microscope cameras often write (near-)linear data; for the images
    tested so far g = 1 put the image colours on the simulated locus far better than 'srgb'."""
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

    alt_value = alt_residual = None
    if reference.system == "MoO3":
        alt_v, alt_d = best_alternative(flat, index, ~bad, reference, ambiguity_gap_nm)
        alt_value, alt_residual = alt_v.reshape(tcd.shape), alt_d.reshape(tcd.shape)

    calibrated_display = np.clip(C.colour.models.eotf_inverse_sRGB(np.clip(lin_cal, 0, 1)), 0, 1)
    settings = dict(max_side=max_side, smooth_sigma=smooth_sigma, max_residual=max_residual,
                    correction_mode=correction_mode, ambiguity_gap_nm=ambiguity_gap_nm, gamma=gamma,
                    image=str(image_path))
    return TCDResult(img, calibrated_display, Jab, tcd, index.reshape(tcd.shape), value.reshape(tcd.shape),
                     residual.reshape(tcd.shape), alt_value, alt_residual, saturated, reference,
                     all_anchors, M, scale, settings)


def layer_classes(result: TCDResult, max_class: int | None = None, clean: int = 5) -> np.ndarray:
    """PS only: integer layer map (0 = substrate, 1 = 1L, ...; -1 unassigned), majority-filtered."""
    cls = np.where(result.assigned, np.floor(np.nan_to_num(result.value, nan=0) + 0.5).astype(int), -1)  # half up
    if max_class is not None:
        cls = np.where(cls > max_class, max_class, cls)
    n = int(np.nanmax(result.reference.value)) + 1
    return mode_filter(cls, clean, n) if clean and clean > 1 else cls


def summary_rows(result: TCDResult) -> list[dict]:
    total = result.index.size
    rows = [dict(item="unassigned", pixels=int((~result.assigned).sum()),
                 fraction=float((~result.assigned).mean()), median_residual=np.nan)]
    if result.reference.system == "PS":
        cls = layer_classes(result, clean=0)
        for c in range(int(np.nanmax(result.reference.value)) + 1):
            m = cls == c
            rows.append(dict(item="substrate" if c == 0 else f"{c}L", pixels=int(m.sum()),
                             fraction=float(m.sum() / total),
                             median_residual=float(np.median(result.residual[m])) if m.any() else np.nan))
    else:
        v = result.value
        flake = result.assigned & (v > 20)
        close = flake & (result.alt_residual - result.residual < 2.0)
        rows.append(dict(item="flake px with alternative thickness within 2 dE", pixels=int(close.sum()),
                         fraction=float(close.sum() / max(flake.sum(), 1)), median_residual=np.nan))
        edges = np.arange(0, np.nanmax(result.reference.value) + 50, 50)
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (v >= lo) & (v < hi)
            rows.append(dict(item=f"{lo:.0f}-{hi:.0f} nm", pixels=int(m.sum()), fraction=float(m.sum() / total),
                             median_residual=float(np.median(result.residual[m])) if m.any() else np.nan))
    return rows
