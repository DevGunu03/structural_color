"""Figures for references and TCD results."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .image_tcd import TCDResult, layer_classes
from .reference import Reference

UNASSIGNED = (0.82, 0.82, 0.82)


def _cbar(fig, im, ax, **kw):
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.05)
    return fig.colorbar(im, cax=cax, **kw)


def reference_figure(ref: Reference, path) -> None:
    """Simulated colour bar and Delta E (CAM02-UCS) from the substrate versus structure."""
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 4.6), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 2.2]))
    x = ref.value
    ax0.imshow(ref.sRGB[None, :, :], aspect="auto", extent=[x.min(), x.max(), 0, 1])
    ax0.set_yticks([])
    ax0.set_title(f"{ref.system}: simulated colour ({ref.meta.get('illuminant')}, NA={ref.meta.get('NA')})",
                  fontsize=10)
    ax1.plot(x, ref.dE_substrate, color="#2a6f97", lw=1.6)
    ax1.set_ylabel("ΔE from substrate (CAM02-UCS)")
    ax1.set_xlabel("MoO$_3$ thickness (nm)" if ref.system == "MoO3" else
                   "effective layer number (integer = dense layer, fraction = packing of top layer)")
    ax1.grid(alpha=0.3)
    fig.text(0.01, 0.005, ref.meta.get("stack", ""), fontsize=7, color="0.4")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
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

    if ref.system == "PS":
        cls = layer_classes(result)
        n = int(np.nanmax(ref.value))
        dense = [ref.find(f"{k}L") if k else 0 for k in range(n + 1)]
        colours = [UNASSIGNED] + [tuple(ref.sRGB[i]) for i in dense]
        cmap = ListedColormap(colours)
        norm = BoundaryNorm(np.arange(-1.5, n + 1.5), cmap.N)
        im = a[3].imshow(cls, cmap=cmap, norm=norm, interpolation="nearest")
        cb = _cbar(fig, im, a[3], ticks=np.arange(-1, n + 1))
        cb.ax.set_yticklabels(["unassigned", "substrate"] + [f"{k}L" for k in range(1, n + 1)])
        a[3].set_title("Layer number (shown in simulated colours)")
    else:
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad(UNASSIGNED)
        im = a[3].imshow(np.ma.masked_invalid(result.value), cmap=cmap, vmin=0, vmax=np.nanmax(ref.value))
        _cbar(fig, im, a[3], label="MoO$_3$ thickness (nm)")
        a[3].set_title("Thickness (nearest reference, grey = unassigned)")

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
