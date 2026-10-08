"""Simulated colour references ("TCD references") built from a system file.

A reference is a table of candidate structures (the sweep of a system file) with the
simulated reflectance colour of each in CIE XYZ, CAM02-UCS, CIELAB and sRGB. Row 0 is the
calibration structure, normally the bare substrate. One column is the *label*, the number the
thickness / layer map reports (thickness in nm, layer number, ...); its description (unit,
title, whether it is rounded to whole classes, the ambiguity gap) is stored in the JSON file
next to the CSV, so the mapping code needs nothing but the reference.

    ref = build_reference(load_system("moo3"))
    ref.save("refs/moo3_D65.csv")          # + refs/moo3_D65.json, readable by the MATLAB code
    ref = Reference.load("refs/moo3_D65.csv")
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import colorimetry as C
from .system import FORMAT_VERSION, System

COLOUR_COLS = ["X", "Y", "Z", "Jp", "ap", "bp", "L", "a", "b", "sR", "sG", "sB", "dE_substrate"]

# label descriptions for references written before system files existed
_LEGACY_LABELS = {
    "MoO3": dict(name="thickness_nm", unit="nm", title="MoO3 thickness", classes=False, class_name="{}",
                 zero_name="", gap=40.0),
    "PS": dict(name="eff_layers", unit="layers", title="Layer number", classes=True, class_name="{}L",
               zero_name="substrate", gap=0.5),
}


@dataclass
class Reference:
    name: str
    params: dict[str, np.ndarray]  # per-row structural parameters, including the label column
    XYZ: np.ndarray
    label: dict
    meta: dict = field(default_factory=dict)
    spectra: np.ndarray | None = None

    def __post_init__(self):
        self.XYZ = np.asarray(self.XYZ, dtype=float)
        self.Jab = C.XYZ_to_cam02ucs(self.XYZ)
        self.Lab = C.XYZ_to_lab(self.XYZ)
        self.sRGB = C.XYZ_to_srgb(self.XYZ)
        self.dE_substrate = C.delta_e(self.Jab, self.Jab[0])

    def __len__(self):
        return len(self.XYZ)

    @property
    def value(self) -> np.ndarray:
        """The label of every row: thickness, layer number, ... (see ``label``)."""
        return self.params[self.label["name"]]

    @property
    def classes(self) -> bool:
        return bool(self.label.get("classes"))

    def axis_title(self) -> str:
        unit, title = self.label.get("unit", ""), self.label.get("title", self.label["name"])
        return title + (f" ({unit})" if unit and unit.lower().rstrip("s") not in title.lower() else "")

    def class_name(self, k: int) -> str:
        if k == 0 and self.label.get("zero_name"):
            return self.label["zero_name"]
        return self.label.get("class_name", "{}").replace("{}", str(k))

    def value_text(self, v: float) -> str:
        unit = self.label.get("unit", "")
        return f"{v:g}" + (f" {unit}" if unit and unit != "layers" else "")

    def class_range(self) -> range:
        """Whole-number classes the label can be rounded to (``classes`` references)."""
        v = self.value
        return range(int(np.floor(np.nanmin(v) + 0.5)), int(np.floor(np.nanmax(v) + 0.5)) + 1)

    def find(self, spec: str) -> int:
        """Row index for an anchor: 'substrate' / 'bare' (row 0), a label value such as '2', '2L',
        '250nm' or '250 nm' (nearest row), or exact parameters 'layers=2,packing=0.5'."""
        s = spec.strip()
        low = s.lower().replace(" ", "")
        if low in ("substrate", "bare", (self.label.get("zero_name") or "substrate").lower()):
            return 0
        if "=" in s:
            mask = np.ones(len(self), bool)
            for part in s.split(","):
                k, v = (x.strip() for x in part.split("=", 1))
                if k not in self.params:
                    raise KeyError(f"'{k}' is not a column of this reference ({', '.join(self.params)})")
                mask &= np.isclose(self.params[k], float(v), atol=1e-6)
            hit = np.flatnonzero(mask)
            if not len(hit):
                raise KeyError(f"No reference row with {spec}")
            return int(hit[0])
        return int(np.argmin(np.abs(self.value - self.parse_value(spec))))

    def parse_value(self, spec: str) -> float:
        """Label value of a text such as '230nm', '230 nm', '2L', '2' or 'substrate'."""
        low = spec.strip().lower().replace(" ", "")
        if low in ("substrate", "bare", (self.label.get("zero_name") or "substrate").lower()):
            return float(self.value[0])
        num = low
        for suffix in (self.label.get("unit", "").lower().replace(" ", ""), "nm", "l"):
            if suffix and num.endswith(suffix):
                num = num[: -len(suffix)]
                break
        try:
            return float(num)
        except ValueError:
            raise KeyError(f"Cannot read '{spec}': use 'substrate', a {self.label['title']} value "
                           f"(e.g. '{self.value_text(self.value[len(self) // 2])}') or 'param=value,...'") from None

    def save(self, path: str | Path, with_spectra: bool = False) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        keys = list(self.params)
        cols = np.column_stack([self.params[k] for k in keys] + [self.XYZ, self.Jab, self.Lab, self.sRGB,
                                                                 self.dE_substrate])
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(keys + COLOUR_COLS)
            for row in cols:
                w.writerow([f"{v:.6g}" for v in row])
        meta = dict(format=FORMAT_VERSION, name=self.name, label=self.label, **self.meta,
                    columns=keys + COLOUR_COLS)
        path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
        if with_spectra and self.spectra is not None:
            header = "wavelength_nm," + ",".join(f"row{i}" for i in range(len(self)))
            np.savetxt(path.with_name(path.stem + "_spectra.csv"), np.column_stack([C.WAVELENGTHS, self.spectra.T]),
                       delimiter=",", fmt="%.6g", header=header, comments="")

    @classmethod
    def load(cls, path: str | Path) -> "Reference":
        path = Path(path)
        meta = json.loads(path.with_suffix(".json").read_text())
        data = np.genfromtxt(path, delimiter=",", names=True)
        names = data.dtype.names
        n_par = len(names) - len(COLOUR_COLS)
        params = {k: np.atleast_1d(np.asarray(data[k], float)) for k in names[:n_par]}
        XYZ = np.column_stack([np.atleast_1d(data[c]) for c in ("X", "Y", "Z")])
        name = meta.pop("name", None) or meta.pop("system")
        meta.pop("system", None)
        meta.pop("columns", None)
        meta.pop("format", None)
        label = meta.pop("label", None) or _LEGACY_LABELS[name]
        return cls(name, params, XYZ, label, meta)

    def subset(self, mask: np.ndarray) -> "Reference":
        mask = np.asarray(mask, bool).copy()
        mask[0] = True  # always keep the calibration row
        spectra = None if self.spectra is None else self.spectra[mask]
        return Reference(self.name, {k: v[mask] for k, v in self.params.items()}, self.XYZ[mask],
                         dict(self.label), dict(self.meta), spectra)


def build_reference(system: System, illuminant: str | None = None, na: float | None = None) -> Reference:
    """Simulate every candidate of ``system``: TMM spectra -> XYZ -> CAM02-UCS, CIELAB, sRGB.

    ``illuminant`` and ``na`` override the system file's optics block.
    """
    illuminant = str(illuminant or system.optics["illuminant"])
    na = system.optics["na"] if na is None else float(na)
    spectra = system.spectra(na=na)
    XYZ = C.reflectance_to_XYZ(spectra, illuminant)
    meta = dict(description=system.description, stack=system.describe(), illuminant=illuminant, NA=na,
                na_weighting=system.optics["na_weighting"], constants=system.constants,
                system_file=system.relative_path(), created_by="python tcd (python/tcd/reference.py)",
                definition=system.definition)
    return Reference(system.name, {k: np.asarray(v, float) for k, v in system.rows.items()}, XYZ,
                     dict(system.label), meta, spectra)


def ambiguity(ref: Reference, gap: float | None = None) -> np.ndarray:
    """For each row, the smallest Delta E (CAM02-UCS) to any row whose label differs by more
    than ``gap``. Small values mean the colour alone cannot tell this structure from another."""
    gap = ref.label["gap"] if gap is None else gap
    out = np.full(len(ref), np.nan)
    for s in range(0, len(ref), 1000):  # chunks keep memory at 1000 x rows
        d = np.linalg.norm(ref.Jab[s:s + 1000, None, :] - ref.Jab[None, :, :], axis=-1)
        d[np.abs(ref.value[s:s + 1000, None] - ref.value[None, :]) <= gap] = np.inf
        m = d.min(axis=1)
        out[s:s + 1000] = np.where(np.isfinite(m), m, np.nan)
    return out
