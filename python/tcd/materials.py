"""Optical constants (n, k) and effective-medium mixing.

The default source is ``data/nk_library.csv``: a CSV export of the
``Index_of_Refraction_library.xls`` used by the MATLAB transfer-matrix code
(TransferMatrix_Updated/), so both codes see identical n, k values.
Two-column or three-column text files (wavelength, n[, k]) are also accepted,
which covers the files in the ``n,k/`` folder (wavelength in um or nm).
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
# repository layout: <repo>/data shared with the MATLAB code; standalone layout: <package parent>/data
DATA_DIR = next((p for p in (_HERE.parents[1] / "data", _HERE.parent / "data") if (p / "nk_library.csv").exists()),
                _HERE.parents[1] / "data")
LIBRARY_CSV = DATA_DIR / "nk_library.csv"

# Close-packed sphere fractions used by the PS-bead model
F_HEX_MONOLAYER = np.pi / (3 * np.sqrt(3))  # 0.6046, one hexagonal layer in a slab of height d
F_CLOSE_PACKED = np.pi / (3 * np.sqrt(2))  # 0.7405, bulk fcc/hcp packing


def _key(name: str) -> str:
    """Normalise a material name so 'Si-Franta', 'Si_Franta' and 'si franta' match."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


@lru_cache(maxsize=None)
def _library(path: str = str(LIBRARY_CSV)) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    raw = np.genfromtxt(path, delimiter=",", names=True)
    cols = raw.dtype.names
    wl = np.asarray(raw[cols[0]], dtype=float)
    out = {}
    for c in cols[1:]:
        if c.endswith("_n"):
            base = c[:-2]
            k_col = base + "_k"
            if k_col in cols:
                out[_key(base)] = (wl, np.asarray(raw[c], float), np.asarray(raw[k_col], float))
    return out


def library_materials() -> list[str]:
    return sorted(_library().keys())


def load_nk_file(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a whitespace/comma separated (wavelength, n[, k]) file.

    Wavelengths below 50 are taken to be in micrometres and converted to nm.
    """
    txt = Path(path).read_text().replace(",", " ")
    data = np.array([[float(v) for v in line.split()] for line in txt.splitlines() if line.strip()
                     and not line.strip()[0].isalpha()])
    wl = data[:, 0] * (1000.0 if data[:, 0].max() < 50 else 1.0)
    n = data[:, 1]
    k = data[:, 2] if data.shape[1] > 2 else np.zeros_like(n)
    order = np.argsort(wl)
    return wl[order], n[order], k[order]


def nk(material: str | complex | float, wavelengths: np.ndarray) -> np.ndarray:
    """Complex refractive index n + ik on ``wavelengths`` (nm).

    ``material`` can be a library name (e.g. 'SiO2-Franta'), a path to an n,k
    text file, or a constant number. Interpolation and extrapolation are linear,
    matching ``interp1(..., 'linear', 'extrap')`` in the MATLAB code.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)
    if isinstance(material, (int, float, complex, np.number)):
        return np.full(wavelengths.shape, complex(material))
    p = Path(str(material))
    if p.suffix and p.exists():
        wl, n, k = load_nk_file(p)
    else:
        lib = _library()
        key = _key(str(material))
        if key not in lib:
            raise KeyError(f"Material '{material}' not in library. Available: {', '.join(library_materials())}")
        wl, n, k = lib[key]
    return _interp_extrap(wl, n, wavelengths) + 1j * _interp_extrap(wl, k, wavelengths)


def _interp_extrap(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    yq = np.interp(xq, x, y)
    lo, hi = xq < x[0], xq > x[-1]
    if lo.any():
        yq[lo] = y[0] + (xq[lo] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    if hi.any():
        yq[hi] = y[-1] + (xq[hi] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return yq


def maxwell_garnett(n_host: np.ndarray, n_incl: np.ndarray, fill: float) -> np.ndarray:
    """Effective index of spherical inclusions (fraction ``fill``) in a host (Maxwell-Garnett)."""
    if not 0.0 <= fill <= 1.0:
        raise ValueError(f"fill fraction {fill} outside [0, 1]")
    eh, ei = n_host**2, n_incl**2
    eps = eh * (2 * (1 - fill) * eh + (1 + 2 * fill) * ei) / ((2 + fill) * eh + (1 - fill) * ei)
    return np.sqrt(eps)
