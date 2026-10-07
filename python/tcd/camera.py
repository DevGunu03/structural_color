"""Measure how a microscope camera encodes intensity, from an exposure series.

Take 4-6 images of the same field (a bare substrate is ideal) with auto-exposure, auto-gain
and auto-white-balance OFF, changing only the exposure time. Recorded intensity is
proportional to exposure time, so the pixel value v follows v ~ t^s. A linear camera gives
s = 1; an sRGB-encoded one gives s ~ 0.42-0.45 over the mid-tones. The decoding exponent for
``run_tcd.py map --gamma`` is g = 1/s.
"""
from __future__ import annotations

import numpy as np

from .image_tcd import Roi, load_image


def region_means(paths: list[str], roi: Roi | None = None, lo: float = 0.02, hi: float = 0.985) -> np.ndarray:
    """Per-channel mean value (0-1) in ``roi`` (default: central half of the frame) for each image;
    NaN where the region is clipped (more than 1 % of its pixels outside [lo, hi])."""
    out = []
    for p in paths:
        img, _ = load_image(p, max_side=None)
        if roi is None:
            h, w = img.shape[:2]
            patch = img[h // 4: 3 * h // 4, w // 4: 3 * w // 4]
        else:
            x, y, rw, rh = roi
            patch = img[y:y + rh, x:x + rw]
        px = patch.reshape(-1, 3)
        clipped = ((px < lo) | (px > hi)).mean(axis=0) > 0.01
        out.append(np.where(clipped, np.nan, px.mean(axis=0)))
    return np.array(out)


def fit_decoding_gamma(values: np.ndarray, exposures: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Log-log slope s per channel (v ~ t^s) and the decoding exponent g = 1/s."""
    t = np.log(np.asarray(exposures, float))
    slopes = []
    for c in range(values.shape[1]):
        v = values[:, c]
        ok = np.isfinite(v)
        slopes.append(np.polyfit(t[ok], np.log(v[ok]), 1)[0] if ok.sum() >= 3 else np.nan)
    s = np.array(slopes)
    return s, 1.0 / s
