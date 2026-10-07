"""Choosing the calibration regions (substrate and optional known layers)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .image_tcd import Roi


def suggest_roi(img: np.ndarray, target_rgb: np.ndarray, size: int = 40, scale: float = 1.0) -> Roi:
    """Square window (in original-image pixels) whose pixels are, on average, closest to
    ``target_rgb`` (sRGB 0-1) - e.g. the substrate colour noted for a given microscope."""
    dist = np.linalg.norm(img - np.asarray(target_rgb)[None, None, :], axis=-1)
    s = max(3, int(round(size * scale)))
    mean = ndimage.uniform_filter(dist, s, mode="nearest")
    half = s // 2
    mean[:half], mean[-half:], mean[:, :half], mean[:, -half:] = np.inf, np.inf, np.inf, np.inf
    y, x = np.unravel_index(np.argmin(mean), mean.shape)
    return (int((x - half) / scale), int((y - half) / scale), size, size)


def dominant_colour(img: np.ndarray, bins: int = 16) -> np.ndarray:
    """Most common colour (coarse 3-D histogram mode); the substrate in sparse-flake images."""
    q = np.clip((img.reshape(-1, 3) * bins).astype(int), 0, bins - 1)
    code = q[:, 0] * bins * bins + q[:, 1] * bins + q[:, 2]
    mode = np.bincount(code).argmax()
    return img.reshape(-1, 3)[code == mode].mean(axis=0)


def select_roi_interactive(img: np.ndarray, title: str, scale: float = 1.0) -> Roi:
    """Click two opposite corners of a region; returns it in original-image pixels."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.imshow(img)
    ax.set_title(title + "\n(click two opposite corners)")
    pts = plt.ginput(2, timeout=0)
    plt.close(fig)
    if len(pts) < 2:
        raise RuntimeError("ROI selection cancelled")
    (x0, y0), (x1, y1) = pts
    x, y = min(x0, x1) / scale, min(y0, y1) / scale
    return (int(x), int(y), max(int(abs(x1 - x0) / scale), 1), max(int(abs(y1 - y0) / scale), 1))
