"""Optical constants (n, k) and effective-medium mixing.

Sources of n + ik (k > 0 means absorption), all returned on the wavelengths asked for (nm):
  * the library ``data/nk_library.csv``, a CSV export of the ``Index_of_Refraction_library.xls``
    used by the original MATLAB transfer-matrix code, so both codes see identical values;
  * a text file of (wavelength, n[, k]) columns: wavelength in nm, or in um if every value
    is below 50; comma, tab or space separated; header lines are skipped;
  * a constant, or a Cauchy or Sellmeier dispersion formula (``cauchy``, ``sellmeier``).

Tabulated data are interpolated and extrapolated linearly, as ``interp1(..., 'linear',
'extrap')`` does in MATLAB. The library covers 300-800 nm, so 800-830 nm is extrapolated.

Mixtures of two materials (beads in air, porous films, rough layers) use an effective-medium
rule: ``maxwell_garnett`` (isolated spheres in a host), ``bruggeman`` (both phases on an equal
footing, e.g. a 50/50 rough interface) or ``linear`` (volume-averaged permittivity).
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
                out[_key(base)] = (wl, np.asarray(raw[c], float), np.asarray(raw[k_col], float), base)
    return out


def library_materials() -> list[str]:
    """Names of the library materials (any spelling that differs only in case, '-', '_' or spaces works)."""
    return sorted((v[3] for v in _library().values()), key=str.lower)


def load_nk_file(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a (wavelength, n[, k]) text file: comma, tab or space separated, header lines skipped.

    Wavelengths are in nm, or in um when every value is below 50 (converted to nm).
    """
    rows = []
    for line in Path(path).read_text().replace(",", " ").replace(";", " ").splitlines():
        try:
            vals = [float(v) for v in line.split()]
        except ValueError:  # header or comment line
            continue
        if len(vals) >= 2:
            rows.append(vals[:3])
    if not rows:
        raise ValueError(f"no (wavelength, n[, k]) rows found in {path}")
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows])
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
        wl, n, k, _ = lib[key]
    return _interp_extrap(wl, n, wavelengths) + 1j * _interp_extrap(wl, k, wavelengths)


def _interp_extrap(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    yq = np.interp(xq, x, y)
    lo, hi = xq < x[0], xq > x[-1]
    if lo.any():
        yq[lo] = y[0] + (xq[lo] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    if hi.any():
        yq[hi] = y[-1] + (xq[hi] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return yq


def cauchy(wavelengths: np.ndarray, A: float, B: float = 0.0, C: float = 0.0) -> np.ndarray:
    """n = A + B / lambda^2 + C / lambda^4 with lambda in micrometres (B in um^2, C in um^4)."""
    lam = np.asarray(wavelengths, dtype=float) / 1000.0
    return A + B / lam**2 + C / lam**4


def sellmeier(wavelengths: np.ndarray, B: list[float], C: list[float]) -> np.ndarray:
    """n^2 = 1 + sum_i B_i lambda^2 / (lambda^2 - C_i), lambda in um, C_i in um^2."""
    if len(B) != len(C):
        raise ValueError("Sellmeier B and C need the same number of terms")
    lam2 = (np.asarray(wavelengths, dtype=float) / 1000.0) ** 2
    eps = 1.0 + sum(b * lam2 / (lam2 - c) for b, c in zip(B, C))
    return np.sqrt(eps + 0j)


def _check_fill(fill: float) -> None:
    if not 0.0 <= fill <= 1.0:
        raise ValueError(f"fill fraction {fill} outside [0, 1]")


def maxwell_garnett(n_host: np.ndarray, n_incl: np.ndarray, fill: float) -> np.ndarray:
    """Effective index of spherical inclusions (volume fraction ``fill``) in a host.

    eps = eps_h (2 (1 - f) eps_h + (1 + 2 f) eps_i) / ((2 + f) eps_h + (1 - f) eps_i),
    the same formula as MaxwellGarnettFormula in TransferMatrix_packing.m. Exact at f = 0
    (host) and f = 1 (inclusion); meant for f well below the percolation of the inclusions.
    """
    _check_fill(fill)
    eh, ei = n_host**2, n_incl**2
    eps = eh * (2 * (1 - fill) * eh + (1 + 2 * fill) * ei) / ((2 + fill) * eh + (1 - fill) * ei)
    return np.sqrt(eps)


def bruggeman(n_host: np.ndarray, n_incl: np.ndarray, fill: float) -> np.ndarray:
    """Bruggeman effective index: f (e_i - e) / (e_i + 2 e) + (1 - f) (e_h - e) / (e_h + 2 e) = 0.

    The quadratic 2 e^2 - b e - e_i e_h = 0, b = (3 f - 1) e_i + (2 - 3 f) e_h, has two roots;
    the physical one has the larger imaginary part (or, for lossless media, the positive one).
    """
    _check_fill(fill)
    eh, ei = n_host**2 + 0j, n_incl**2 + 0j
    b = (3 * fill - 1) * ei + (2 - 3 * fill) * eh
    root = np.sqrt(b**2 + 8 * ei * eh)
    e1, e2 = (b + root) / 4, (b - root) / 4
    pick1 = (e1.imag > e2.imag + 1e-12) | ((np.abs(e1.imag - e2.imag) <= 1e-12) & (e1.real >= e2.real))
    return np.sqrt(np.where(pick1, e1, e2))


def linear_mix(n_host: np.ndarray, n_incl: np.ndarray, fill: float) -> np.ndarray:
    """Volume-averaged permittivity: eps = f eps_i + (1 - f) eps_h."""
    _check_fill(fill)
    return np.sqrt(fill * n_incl**2 + (1 - fill) * n_host**2 + 0j)


EMA_RULES = {"maxwell-garnett": maxwell_garnett, "bruggeman": bruggeman, "linear": linear_mix}
