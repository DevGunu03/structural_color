"""Optical constants (n, k) and effective-medium mixing.

Sources of n + ik (k > 0 means absorption), all returned on the wavelengths asked for (nm):
  * the library: ``data/nk_library.csv`` (a CSV export of the ``Index_of_Refraction_library.xls``
    used by the original MATLAB transfer-matrix code, so both codes see identical values) plus
    one file per added material in ``data/materials/`` (``add_material`` / ``add-material``);
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
from dataclasses import dataclass, field
from datetime import date
from functools import lru_cache
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
# repository layout: <repo>/data shared with the MATLAB code; standalone layout: <package parent>/data
DATA_DIR = next((p for p in (_HERE.parents[1] / "data", _HERE.parent / "data") if (p / "nk_library.csv").exists()),
                _HERE.parents[1] / "data")
LIBRARY_CSV = DATA_DIR / "nk_library.csv"
MATERIALS_DIR = DATA_DIR / "materials"            # added materials: one file per material, part of the library
NK_SUFFIXES = (".csv", ".txt", ".dat", ".nk", ".tsv")
_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")

# Close-packed sphere fractions used by the PS-bead model
F_HEX_MONOLAYER = np.pi / (3 * np.sqrt(3))  # 0.6046, one hexagonal layer in a slab of height d
F_CLOSE_PACKED = np.pi / (3 * np.sqrt(2))  # 0.7405, bulk fcc/hcp packing


@dataclass
class Material:
    name: str
    wl: np.ndarray          # nm, increasing
    n: np.ndarray
    k: np.ndarray
    origin: str             # 'nk_library.csv' or the file in data/materials/
    meta: dict = field(default_factory=dict)   # source, reference, notes, added (from '# key: value' lines)

    @property
    def source(self) -> str:
        return self.meta.get("source", "")


def _key(name: str) -> str:
    """Normalise a material name so 'Si-Franta', 'Si_Franta' and 'si franta' match."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def parse_nk_text(text: str, units: str = "auto", where: str = "file") -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """(wavelength nm, n, k, metadata) from the text of an n,k file.

    Understands:
      * columns wavelength, n[, k] separated by commas, tabs, spaces or semicolons;
      * header lines (any line that is not numbers) and '# key: value' metadata lines;
      * the refractiveindex.info CSV export, where an 'wl,n' block is followed by an 'wl,k' block.
    Wavelengths are nm, or um when every value is below 50 (``units='auto'``); 'nm' or 'um'
    forces the unit.
    """
    meta, blocks, cur = {}, [], None
    for line in text.splitlines():
        s = line.strip().lstrip("﻿")
        if not s:
            continue
        if s.startswith("#"):
            m = re.match(r"#\s*([A-Za-z][A-Za-z _]*?)\s*:\s*(.*)$", s)
            if m:
                meta[m.group(1).strip().lower()] = m.group(2).strip()
            continue
        toks = [t for t in re.split(r"[,\s;]+", s) if t]
        try:
            vals = [float(t) for t in toks]
        except ValueError:          # a header line starts a new block
            cur = dict(header=s.lower(), rows=[])
            blocks.append(cur)
            continue
        if cur is None:
            cur = dict(header="", rows=[])
            blocks.append(cur)
        cur["rows"].append(vals)
    blocks = [b for b in blocks if b["rows"]]
    if not blocks:
        raise ValueError(f"{where}: no (wavelength, n[, k]) rows found")

    def table(b):
        w = min(len(r) for r in b["rows"])
        if w < 2:
            raise ValueError(f"{where}: rows need at least two columns (wavelength, n)")
        return np.array([r[:w] for r in b["rows"]], dtype=float)

    if len(blocks) == 1:
        d = table(blocks[0])
        wl, n = d[:, 0], d[:, 1]
        k = d[:, 2] if d.shape[1] > 2 else np.zeros_like(n)
    elif len(blocks) == 2:
        d1, d2 = table(blocks[0]), table(blocks[1])
        is_k = [bool(re.search(r"(^|[^a-z])k([^a-z]|$)", b["header"])) for b in blocks]
        dn, dk = (d2, d1) if is_k[0] and not is_k[1] else (d1, d2)
        wl, n = dn[:, 0], dn[:, 1]
        o = np.argsort(dk[:, 0])
        k = np.interp(wl, dk[o, 0], dk[o, 1])
    else:
        raise ValueError(f"{where}: found {len(blocks)} separate tables; expected one (wavelength, n, k) table "
                         f"or an n table followed by a k table")
    if units == "um" or (units == "auto" and wl.max() < 50):
        wl = wl * 1000.0
    elif units not in ("auto", "nm"):
        raise ValueError("units must be 'auto', 'nm' or 'um'")
    o = np.argsort(wl)
    wl, n, k = wl[o], n[o], k[o]
    if np.any(np.diff(wl) == 0):
        raise ValueError(f"{where}: the same wavelength appears twice")
    return wl, n, k, meta


def read_nk_file(path: str | Path, units: str = "auto") -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    return parse_nk_text(Path(path).read_text(encoding="utf-8", errors="replace"), units, str(path))


def load_nk_file(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(wavelength nm, n, k) of an n,k file; see ``parse_nk_text`` for the formats."""
    wl, n, k, _ = read_nk_file(path)
    return wl, n, k


@lru_cache(maxsize=None)
def _library() -> dict[str, Material]:
    """Library materials by normalised name: nk_library.csv first, then the files in data/materials/."""
    raw = np.genfromtxt(LIBRARY_CSV, delimiter=",", names=True)
    cols = raw.dtype.names
    wl = np.asarray(raw[cols[0]], dtype=float)
    out = {}
    for c in cols[1:]:
        if c.endswith("_n") and c[:-2] + "_k" in cols:
            base = c[:-2]
            out[_key(base)] = Material(base, wl, np.asarray(raw[c], float), np.asarray(raw[base + "_k"], float),
                                       "data/nk_library.csv", {"source": "nk_library.csv (Index_of_Refraction_library.xls)"})
    if MATERIALS_DIR.is_dir():
        for f in sorted(MATERIALS_DIR.iterdir()):
            if f.suffix.lower() not in NK_SUFFIXES or _key(f.stem) in out:
                continue                      # a built-in name always wins; add_material refuses clashes
            try:
                w, n, k, meta = read_nk_file(f)
            except ValueError:
                continue
            out[_key(f.stem)] = Material(meta.get("name", f.stem), w, n, k, f"data/materials/{f.name}", meta)
    return out


def library_materials() -> list[str]:
    """Names of the library materials (any spelling that differs only in case, '-', '_' or spaces works)."""
    return sorted((m.name for m in _library().values()), key=str.lower)


def library_entries() -> list[Material]:
    return sorted(_library().values(), key=lambda m: m.name.lower())


def add_material(name: str, source_file: str | Path, source: str, reference: str = "", notes: str = "",
                 units: str = "auto", replace: bool = False, folder: Path | None = None) -> tuple[Path, list[str]]:
    """Check an n,k file and save it as ``data/materials/<name>.csv``, part of the library from then on.

    Returns the new file and a list of warnings (e.g. wavelengths that will be extrapolated).
    Refuses a name that the library already has, unless ``replace`` and the existing entry is an
    added material (built-in entries of nk_library.csv are never replaced).
    """
    folder = Path(folder) if folder else MATERIALS_DIR
    if not _NAME.match(name):
        raise ValueError(f"name '{name}': use letters, digits, '-' and '_', starting with a letter (e.g. MoS2_bulk)")
    if not source.strip():
        raise ValueError("give the source of the data (paper, database or measurement), so others can trust it")
    existing = _library().get(_key(name))
    target = folder / f"{name}.csv"
    if existing is not None:
        if existing.origin == "data/nk_library.csv":
            raise ValueError(f"'{name}' is already a built-in library material ({existing.name}); choose another name")
        if not replace:
            raise ValueError(f"'{name}' already exists ({existing.origin}); add --replace to overwrite it")
    wl, n, k, _ = read_nk_file(source_file, units)
    problems = []
    if len(wl) < 2:
        problems.append("need at least two wavelengths")
    if not np.all(np.isfinite(np.r_[wl, n, k])):
        problems.append("the table contains empty or non-numeric values")
    if np.any(n <= 0):
        problems.append("n must be positive")
    if np.any(k < 0):
        problems.append("k must be zero or positive (this code uses n + ik, with k > 0 for absorption)")
    if problems:
        raise ValueError(f"{source_file}: " + "; ".join(problems))
    warnings = []
    if wl[0] > 360:
        warnings.append(f"data start at {wl[0]:g} nm: 360-{wl[0]:g} nm will be extrapolated")
    if wl[-1] < 830:
        warnings.append(f"data end at {wl[-1]:g} nm: {wl[-1]:g}-830 nm will be extrapolated")
    if wl[0] > 400 or wl[-1] < 700:
        warnings.append("the data do not cover 400-700 nm: the simulated colours will be unreliable")
    folder.mkdir(parents=True, exist_ok=True)
    lines = [f"# name: {name}", f"# source: {source.strip()}"]
    if reference.strip():
        lines.append(f"# reference: {reference.strip()}")
    if notes.strip():
        lines.append(f"# notes: {notes.strip()}")
    lines += [f"# added: {date.today().isoformat()}", f"# imported from: {Path(source_file).name}",
              "wavelength_nm,n,k"]
    lines += [f"{a:.6g},{b:.6g},{c:.6g}" for a, b, c in zip(wl, n, k)]
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _library.cache_clear()
    return target, warnings


def nk(material: str | complex | float, wavelengths: np.ndarray) -> np.ndarray:
    """Complex refractive index n + ik on ``wavelengths`` (nm).

    ``material`` can be a library name (e.g. 'SiO2-Franta', or an added material such as
    'MoS2'), a path to an n,k text file, or a constant number. Interpolation and extrapolation
    are linear, matching ``interp1(..., 'linear', 'extrap')`` in the MATLAB code.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)
    if isinstance(material, (int, float, complex, np.number)):
        return np.full(wavelengths.shape, complex(material))
    p = Path(str(material))
    if p.suffix.lower() in NK_SUFFIXES and p.exists():
        wl, n, k = load_nk_file(p)
    else:
        m = _library().get(_key(str(material)))
        if m is None:
            raise KeyError(f"Material '{material}' not in library. Available: {', '.join(library_materials())}")
        wl, n, k = m.wl, m.n, m.k
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
