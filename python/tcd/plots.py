"""Figures: the reference sheet of a system, spectra of chosen structures, and TCD results."""
from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .colorimetry import WAVELENGTHS
from .image_tcd import TCDResult, layer_classes, nice_step
from .reference import Reference, ambiguity

UNASSIGNED = (0.82, 0.82, 0.82)


def _cbar(fig, im, ax, **kw):
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.05)
    return fig.colorbar(im, cax=cax, **kw)


def _ink(rgb) -> str:
    """Black or white text, whichever reads better on ``rgb``."""
    return "black" if np.dot(rgb, [0.299, 0.587, 0.114]) > 0.5 else "white"


def chart_rows(ref: Reference, n: int = 13) -> list[int]:
    """Rows shown as swatches: every class, or about ``n`` evenly spaced label values."""
    if ref.classes:
        return [ref.find(str(c)) for c in ref.class_range()]
    v = ref.value
    w = nice_step(float(np.nanmax(v) - np.nanmin(v)), n - 1)
    targets = np.arange(np.ceil(np.nanmin(v) / w) * w, np.nanmax(v) + w * 1e-6, w)
    return list(dict.fromkeys(int(np.argmin(np.abs(v - t))) for t in targets))


def _with_breaks(x: np.ndarray, *ys: np.ndarray):
    """Insert NaN where the sorted label jumps (a sweep gap) so lines are not drawn across it."""
    steps = np.diff(x)
    typical = np.median(steps[steps > 0]) if np.any(steps > 0) else 1.0
    cut = np.flatnonzero(steps > 5 * typical) + 1
    return [np.insert(np.asarray(v, float), cut, np.nan) for v in (x, *ys)]


def reference_figure(ref: Reference, path) -> None:
    """Reference sheet: what the colours look like, how they change, and where they repeat."""
    order = np.argsort(ref.value, kind="stable")
    x = ref.value[order]
    amb = ambiguity(ref)
    steps = np.diff(x)
    typical = np.median(steps[steps > 0]) if np.any(steps > 0) else 1.0
    fig = plt.figure(figsize=(12, 9.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.7, 2.0, 2.6, 0.25], width_ratios=[1, 1.15],
                          hspace=0.42, wspace=0.22)

    ax0 = fig.add_subplot(gs[0, :])
    xs = np.linspace(x.min(), x.max(), 1500)  # resample on the value axis: rows need not be evenly spaced
    j = np.clip(np.searchsorted(x, xs), 1, len(x) - 1) if len(x) > 1 else np.zeros(len(xs), int)
    if len(x) > 1:
        j = np.where(np.abs(x[j - 1] - xs) <= np.abs(x[j] - xs), j - 1, j)
    strip = ref.sRGB[order][j].copy()
    strip[np.abs(x[j] - xs) > 2.5 * typical] = 1.0  # white where the sweep has no candidates
    ax0.imshow(strip[None, :, :], aspect="auto", extent=[x.min(), x.max(), 0, 1], interpolation="nearest")
    ax0.set_yticks([])
    ax0.set_title(f"{ref.name}: simulated colour under {ref.meta.get('illuminant', 'D65')}, "
                  f"NA = {ref.meta.get('NA', 0):g}" + (f"  -  {ref.meta['description']}"
                                                         if ref.meta.get("description") else ""), fontsize=10)

    ax1 = fig.add_subplot(gs[1, :], sharex=ax0)
    xb, de_b, amb_b, a_b, b_b = _with_breaks(x, ref.dE_substrate[order], amb[order], ref.Jab[order, 1],
                                             ref.Jab[order, 2])
    ax1.plot(xb, de_b, color="#2a6f97", lw=1.6, label="ΔE from the calibration row (row 0)")
    ax1.plot(xb, amb_b, color="#c44e52", lw=1.2,
             label=f"ΔE to the nearest look-alike more than {ref.label['gap']:g} {ref.label.get('unit', '')} away")
    ax1.axhspan(0, 2, color="#c44e52", alpha=0.10, lw=0, label="below ~2 ΔE: hard to tell apart in a micrograph")
    ax1.set_ylabel("ΔE (CAM02-UCS)")
    ax1.set_xlabel(ref.axis_title())
    ax1.set_ylim(bottom=0)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8, loc="best")

    ax2 = fig.add_subplot(gs[2, 0])
    sc = ax2.scatter(ref.Jab[order, 1], ref.Jab[order, 2], c=x, cmap="viridis", s=7, zorder=3)
    ax2.plot(a_b, b_b, color="0.75", lw=0.6, zorder=2)
    ax2.plot(ref.Jab[0, 1], ref.Jab[0, 2], "o", mfc="none", mec="red", ms=10, zorder=4, label="row 0")
    _cbar(fig, sc, ax2, label=ref.axis_title())
    ax2.set_xlabel("a′")
    ax2.set_ylabel("b′")
    ax2.set_aspect("equal", adjustable="datalim")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc="best")
    ax2.set_title("Colour path in CAM02-UCS (loops = repeating colours)", fontsize=10)

    ax3 = fig.add_subplot(gs[2, 1])
    rows = chart_rows(ref)
    ncol = 4 if len(rows) > 8 else max(len(rows), 1)
    nrow = int(np.ceil(len(rows) / ncol))
    for j, r in enumerate(rows):
        cx, cy = j % ncol, nrow - 1 - j // ncol
        ax3.add_patch(Rectangle((cx + 0.04, cy + 0.06), 0.92, 0.88, color=ref.sRGB[r]))
        name = ref.class_name(int(round(ref.value[r]))) if ref.classes else ref.value_text(ref.value[r])
        ax3.text(cx + 0.5, cy + 0.5, name, ha="center", va="center", fontsize=9, color=_ink(ref.sRGB[r]))
    ax3.set_xlim(0, ncol)
    ax3.set_ylim(0, nrow)
    ax3.set_aspect("equal")
    ax3.axis("off")
    ax3.set_title("Colour chart (sRGB, D65 display)", fontsize=10)

    ax4 = fig.add_subplot(gs[3, :])
    ax4.axis("off")
    ax4.text(0, 1, "\n".join(textwrap.wrap(ref.meta.get("stack", ""), 170)), fontsize=7, color="0.35",
             va="top", family="monospace")
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def spectra_figure(spectra: np.ndarray, sRGB: np.ndarray, names: list[str], title: str, path) -> None:
    """Reflectance spectra of a few structures with their simulated colours."""
    fig, (ax, axs) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw=dict(width_ratios=[3, 1]))
    for R, rgb, name in zip(spectra, sRGB, names):
        ax.plot(WAVELENGTHS, R, lw=1.6, color=np.clip(rgb * 0.85, 0, 1), label=name)
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("reflectance")
    ax.set_xlim(WAVELENGTHS[0], WAVELENGTHS[-1])
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    for j, (rgb, name) in enumerate(zip(sRGB, names)):
        y = len(names) - 1 - j
        axs.add_patch(Rectangle((0, y + 0.08), 1, 0.84, color=rgb))
        axs.text(0.5, y + 0.5, name, ha="center", va="center", fontsize=8, color=_ink(rgb))
    axs.set_xlim(0, 1)
    axs.set_ylim(0, len(names))
    axs.axis("off")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _draw_rois(ax, result: TCDResult):
    for a in result.anchors:
        x, y, w, h = (v * result.scale for v in a.roi)
        ax.add_patch(Rectangle((x, y), w, h, fill=False, ec="white", lw=1.5))
        ax.add_patch(Rectangle((x, y), w, h, fill=False, ec="black", lw=0.6, ls="--"))
        ax.text(x, y - 3, a.name, color="white", fontsize=7, va="bottom",
                bbox=dict(fc="black", alpha=0.5, pad=1, lw=0))


def result_figure(result: TCDResult, path, title: str = "") -> None:
    ref = result.reference
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
    a = axes.ravel()

    a[0].imshow(result.image)
    _draw_rois(a[0], result)
    a[0].set_title("Micrograph + calibration regions")
    a[1].imshow(result.calibrated)
    a[1].set_title("Calibrated to simulated substrate colour")

    vmax = np.nanpercentile(result.tcd, 99.5)
    im = a[2].imshow(result.tcd, cmap="cividis", vmin=0, vmax=vmax)
    _cbar(fig, im, a[2], label="ΔE from substrate (CAM02-UCS)")
    a[2].set_title("TCD map (absolute ΔE)")

    if ref.classes:
        cls = layer_classes(result)
        ks = list(ref.class_range())
        colours = [UNASSIGNED] + [tuple(ref.sRGB[ref.find(str(k))]) for k in ks]
        cmap = ListedColormap(colours)
        norm = BoundaryNorm(np.arange(ks[0] - 1.5, ks[-1] + 1.5), cmap.N)
        im = a[3].imshow(cls, cmap=cmap, norm=norm, interpolation="nearest")
        cb = _cbar(fig, im, a[3], ticks=np.arange(ks[0] - 1, ks[-1] + 1))
        cb.ax.set_yticklabels(["unassigned"] + [ref.class_name(k) for k in ks])
        a[3].set_title(f"{ref.label.get('title', 'Class')} (shown in simulated colours)")
    else:
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad(UNASSIGNED)
        im = a[3].imshow(np.ma.masked_invalid(result.value), cmap=cmap, vmin=np.nanmin(ref.value),
                         vmax=np.nanmax(ref.value))
        _cbar(fig, im, a[3], label=ref.axis_title())
        a[3].set_title(f"{ref.label.get('title', 'Value')} (nearest reference, grey = unassigned)")

    cmap = plt.get_cmap("magma").copy()
    im = a[4].imshow(result.residual, cmap=cmap, vmin=0, vmax=max(result.settings["max_residual"] * 1.5, 1))
    _cbar(fig, im, a[4], label="distance to nearest reference (ΔE)")
    a[4].set_title(f"Match residual (unassigned above {result.settings['max_residual']:g})")

    ax = a[5]
    pts = result.Jab.reshape(-1, 3)
    step = max(1, len(pts) // 60000)
    ax.hexbin(pts[::step, 1], pts[::step, 2], gridsize=90, bins="log", cmap="Greys", mincnt=1)
    ax.scatter(ref.Jab[:, 1], ref.Jab[:, 2], c=ref.value, cmap="viridis", s=6, zorder=3)
    ax.plot(ref.Jab[0, 1], ref.Jab[0, 2], "o", mfc="none", mec="red", ms=10, zorder=4, label="substrate")
    for an in result.anchors[1:]:
        ax.plot(*ref.Jab[an.ref_index, 1:], "s", mfc="none", mec="orange", ms=9, zorder=4)
    ax.set_xlabel("a′")
    ax.set_ylabel("b′")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Image colours (grey) vs simulated locus (coloured)")

    for axx in a[:5]:
        axx.set_xticks([])
        axx.set_yticks([])
    fig.suptitle(title or result.settings.get("image", ""), fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
