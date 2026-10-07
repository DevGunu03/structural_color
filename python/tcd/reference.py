"""Simulated colour references ("TCD references") for MoO3 flakes and PS-bead layers.

Each reference is a table of candidate structures (a MoO3 thickness, or a PS layer
count with the packing of its top layer) together with their simulated reflectance
colour in CAM02-UCS, CIELAB and sRGB. Row 0 is always the bare substrate.

Stacks (light incident from the left):
    MoO3 : Air | MoO3 (t)                         | SiO2 (oxide) | Si
    PS   : Air | PS layer N (top, packing p) ... 1 | SiO2 (oxide) | Si
PS layers are Maxwell-Garnett slabs of beads in air.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import colorimetry as C
from .materials import F_CLOSE_PACKED, F_HEX_MONOLAYER, maxwell_garnett, nk
from .tmm import reflectance_na

COLOUR_COLS = ["X", "Y", "Z", "Jp", "ap", "bp", "L", "a", "b", "sR", "sG", "sB", "dE_substrate"]


@dataclass
class Reference:
    system: str
    params: dict[str, np.ndarray]  # per-row structural parameters
    XYZ: np.ndarray
    meta: dict = field(default_factory=dict)
    spectra: np.ndarray | None = None

    def __post_init__(self):
        self.Jab = C.XYZ_to_cam02ucs(self.XYZ)
        self.Lab = C.XYZ_to_lab(self.XYZ)
        self.sRGB = C.XYZ_to_srgb(self.XYZ)
        self.dE_substrate = C.delta_e(self.Jab, self.Jab[0])

    def __len__(self):
        return len(self.XYZ)

    @property
    def value(self) -> np.ndarray:
        """Scalar structural coordinate: thickness (nm) for MoO3, effective layer number for PS."""
        return self.params["thickness_nm"] if self.system == "MoO3" else self.params["eff_layers"]

    def labels(self) -> list[str]:
        if self.system == "MoO3":
            return [f"{t:g} nm" for t in self.params["thickness_nm"]]
        return ["substrate" if n == 0 else f"{int(n)}L ({p:.0%})"
                for n, p in zip(self.params["layers"], self.params["packing"])]

    def find(self, label: str) -> int:
        """Row index for an anchor label: 'substrate', '1L', '2L', ... (dense layer) or '250nm'."""
        s = label.strip().lower().replace(" ", "")
        if s in ("substrate", "0", "0l", "0nm"):
            return 0
        if self.system == "PS" and s.endswith("l"):
            n = int(s[:-1])
            hit = np.flatnonzero((self.params["layers"] == n) & np.isclose(self.params["packing"], 1.0))
        else:
            t = float(s.removesuffix("nm"))
            hit = [int(np.argmin(np.abs(self.value - t)))]
        if len(hit) == 0:
            raise KeyError(f"No reference row for '{label}'")
        return int(hit[0])

    def save(self, path: str | Path, with_spectra: bool = False) -> None:
        path = Path(path)
        keys = list(self.params)
        cols = np.column_stack([self.params[k] for k in keys] + [self.XYZ, self.Jab, self.Lab, self.sRGB,
                                                                 self.dE_substrate])
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(keys + COLOUR_COLS)
            for row in cols:
                w.writerow([f"{v:.6g}" for v in row])
        meta = dict(self.meta, system=self.system, columns=keys + COLOUR_COLS)
        path.with_suffix(".json").write_text(json.dumps(meta, indent=2))
        if with_spectra and self.spectra is not None:
            np.savetxt(path.with_name(path.stem + "_spectra.csv"),
                       np.column_stack([C.WAVELENGTHS, self.spectra.T]), delimiter=",", fmt="%.6g")

    @classmethod
    def load(cls, path: str | Path) -> "Reference":
        path = Path(path)
        meta = json.loads(path.with_suffix(".json").read_text())
        data = np.genfromtxt(path, delimiter=",", names=True)
        names = data.dtype.names
        n_par = len(names) - len(COLOUR_COLS)
        params = {k: np.asarray(data[k], float) for k in names[:n_par]}
        XYZ = np.column_stack([data["X"], data["Y"], data["Z"]])
        return cls(meta.pop("system"), params, XYZ, meta)

    def subset(self, mask: np.ndarray) -> "Reference":
        mask = np.asarray(mask, bool).copy()
        mask[0] = True  # always keep the substrate
        spectra = None if self.spectra is None else self.spectra[mask]
        return Reference(self.system, {k: v[mask] for k, v in self.params.items()}, self.XYZ[mask],
                         dict(self.meta), spectra)


def _colours(spectra, illuminant):
    return C.reflectance_to_XYZ(spectra, illuminant)


def build_moo3(t_max: float = 600, step: float = 1, oxide_nm: float = 100, moo3: str = "aMoO3",
               sio2: str = "SiO2-Franta", si: str = "Si-Franta", illuminant: str = "D65",
               na: float = 0.0) -> Reference:
    """MoO3 thickness series on SiO2/Si. 'aMoO3' is the library's alpha-MoO3 (Lajaunie) data."""
    wl = C.WAVELENGTHS
    n_air, n_moo3, n_ox, n_si = nk("Air", wl), nk(moo3, wl), nk(sio2, wl), nk(si, wl)
    thick = np.arange(0, t_max + step / 2, step, dtype=float)
    spectra = np.array([reflectance_na([n_air, n_moo3, n_ox, n_si], [0, t, oxide_nm, 0], wl, na)
                        for t in thick])
    meta = dict(stack=f"Air | {moo3} (t) | {sio2} {oxide_nm} nm | {si}", illuminant=str(illuminant),
                NA=na, t_max=t_max, step=step)
    return Reference("MoO3", {"thickness_nm": thick}, _colours(spectra, illuminant), meta, spectra)


def ps_stack(n_layers: int, packing: float, bead_nm: float, model: str):
    """(fill fractions, thicknesses) of the PS slabs, top layer first."""
    if n_layers == 0:
        return [], []
    fills, thick = [], []
    for k in range(n_layers, 0, -1):  # k = layer index from the substrate (1 = bottom)
        if model == "slab":  # legacy: every layer a d-thick slab at monolayer packing
            f, d = F_HEX_MONOLAYER, bead_nm
        elif model == "hcp":  # close-packed stacking: layer pitch d*sqrt(2/3)
            f, d = (F_HEX_MONOLAYER, bead_nm) if k == 1 else (F_CLOSE_PACKED, bead_nm * np.sqrt(2 / 3))
        else:
            raise ValueError("model must be 'slab' or 'hcp'")
        fills.append(f * (packing if k == n_layers else 1.0))
        thick.append(d)
    return fills, thick


def build_ps(max_layers: int = 5, packing_step: float = 0.01, bead_nm: float = 300, oxide_nm: float = 100,
             model: str = "slab", ps: str = "PS-beads", sio2: str = "SiO2-Franta", si: str = "Si-Franta",
             illuminant: str = "D65", na: float = 0.0) -> Reference:
    """PS bead layers: bare substrate, then for N = 1..max_layers a top layer whose packing
    goes from ~0 to 1 (fraction of a close-packed layer) on N-1 dense layers."""
    wl = C.WAVELENGTHS
    n_air, n_ps, n_ox, n_si = nk("Air", wl), nk(ps, wl), nk(sio2, wl), nk(si, wl)
    packs = np.round(np.arange(packing_step, 1 + packing_step / 2, packing_step), 6)
    rows = [(0, 0.0)] + [(n, p) for n in range(1, max_layers + 1) for p in packs]
    spectra = []
    for n, p in rows:
        fills, thick = ps_stack(n, p, bead_nm, model)
        layers = [n_air] + [maxwell_garnett(n_air, n_ps, f) for f in fills] + [n_ox, n_si]
        spectra.append(reflectance_na(layers, [0] + thick + [oxide_nm, 0], wl, na))
    layers_col = np.array([r[0] for r in rows], float)
    pack_col = np.array([r[1] for r in rows], float)
    params = {"layers": layers_col, "packing": pack_col,
              "eff_layers": np.where(layers_col == 0, 0, layers_col - 1 + pack_col)}
    meta = dict(stack=f"Air | {ps} layers (d={bead_nm} nm, model={model}) | {sio2} {oxide_nm} nm | {si}",
                illuminant=str(illuminant), NA=na, max_layers=max_layers, packing_step=packing_step)
    return Reference("PS", params, _colours(np.array(spectra), illuminant), meta, np.array(spectra))
