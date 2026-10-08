"""System files: describe a sample once, get its simulated colour reference.

A system file (``systems/*.jsonc``, JSON with ``//`` comments) says what the sample is made of
and which structural parameters vary from place to place. Every combination of the swept
parameters is one candidate structure; the transfer-matrix code turns each into a reflectance
spectrum and the colorimetry code into a colour. The table of candidates and colours is the
reference that ``image_tcd.analyse`` matches micrograph pixels against.

    {
      "name": "moo3",
      "constants": {"oxide_nm": 100},                          // fixed numbers, overridable
      "sweep": {"thickness_nm": {"from": 0, "to": 600, "step": 1}},
      "label": {"name": "thickness_nm", "unit": "nm"},         // what the map reports
      "ambient": "Air",
      "layers": [                                             // light side first
        {"name": "flake", "material": "aMoO3", "thickness": "thickness_nm"},
        {"name": "oxide", "material": "SiO2-Franta", "thickness": "oxide_nm"}
      ],
      "substrate": "Si-Franta"
    }

``systems/README.md`` documents every key; ``systems/template.jsonc`` is a commented starting
point. The MATLAB twin is ``tcd_load_system.m`` and reads the same files.
"""
from __future__ import annotations

import copy
import json
import re
from pathlib import Path

import numpy as np

from . import materials as M
from .expr import ExpressionError, evaluate, names_in
from .tmm import reflectance_na

REPO = Path(__file__).resolve().parents[2]
SYSTEMS_DIR = REPO / "systems"
FORMAT_VERSION = 1

_TOP = {"name", "description", "constants", "materials", "ambient", "layers", "substrate", "sweep", "label",
        "optics", "notes"}
_LAYER = {"name", "material", "thickness", "repeat", "notes"}
_GROUP = {"name", "layers", "repeat", "notes"}
_LABEL = {"name", "value", "unit", "title", "classes", "class_name", "zero_name", "gap", "notes"}
_OPTICS = {"illuminant", "na", "na_weighting", "angles", "notes"}
_MATERIAL_KINDS = {
    "library": {"library"},
    "file": {"file"},
    "n": {"n", "k"},
    "cauchy": {"cauchy", "k"},
    "sellmeier": {"sellmeier", "k"},
    "ema": {"ema", "host", "inclusion", "fraction"},
}
_NAME = re.compile(r"^[A-Za-z_]\w*$")
_RESERVED = {"pi", "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "asin", "acos", "atan", "abs", "floor",
             "ceil", "round", "min", "max"}


class SystemError(ValueError):
    """A problem in a system file, with the location in the file in the message."""


# --- reading ---------------------------------------------------------------------------------
def strip_comments(text: str) -> str:
    """Remove ``//`` comments (outside strings) so the rest parses as plain JSON."""
    out, i, n, in_str = [], 0, len(text), False
    while i < n:
        c = text[i]
        if in_str:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 1
            elif c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            out.append(c)
        elif c == "/" and text[i + 1:i + 2] == "/":
            while i < n and text[i] != "\n":
                i += 1
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def read_system_file(path: str | Path) -> dict:
    text = Path(path).read_text(encoding="utf-8")
    try:
        return json.loads(strip_comments(text))
    except json.JSONDecodeError as e:
        raise SystemError(f"{path}: not valid JSON after removing // comments: {e}") from None


def find_system_file(name: str | Path) -> Path:
    """A path, or the name of a bundled system ('moo3' -> systems/moo3.jsonc)."""
    p = Path(name)
    if p.suffix and p.exists():
        return p
    for cand in (SYSTEMS_DIR / f"{name}.jsonc", SYSTEMS_DIR / f"{name}.json", SYSTEMS_DIR / f"{str(name).lower()}.jsonc"):
        if cand.is_file():
            return cand
    known = ", ".join(sorted(f.stem for f in SYSTEMS_DIR.glob("*.json*")))
    raise FileNotFoundError(f"No system file '{name}' (bundled systems: {known})")


def parse_overrides(items: list[str] | None) -> dict[str, str]:
    """['oxide_nm=285', 'bead_nm=500'] -> {'oxide_nm': '285', 'bead_nm': '500'}."""
    out = {}
    for item in items or []:
        if "=" not in item:
            raise SystemError(f"--set expects name=value, got '{item}'")
        k, v = item.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def load_system(source: str | Path | dict, overrides: dict | None = None) -> "System":
    """Load a system from a file path, a bundled name ('moo3', 'ps', ...) or a dict."""
    if isinstance(source, dict):
        return System(copy.deepcopy(source), None, overrides)
    path = find_system_file(source)
    return System(read_system_file(path), path, overrides)


# --- the system ------------------------------------------------------------------------------
class System:
    """A parsed and checked system file, with its candidate structures expanded into ``rows``."""

    def __init__(self, definition: dict, path: Path | None = None, overrides: dict | None = None):
        self.definition = definition
        self.path = Path(path) if path else None
        self.where = str(path) if path else "system"
        d = definition
        if not isinstance(d, dict):
            raise SystemError(f"{self.where}: the file must contain one JSON object {{...}}")
        self._check_keys(d, _TOP, "top level")
        for key in ("name", "layers", "substrate", "sweep"):
            if key not in d:
                raise SystemError(f"{self.where}: missing required key '{key}'")
        self.name = str(d["name"])
        self.description = str(d.get("description", ""))
        self.constants = self._constants(d.get("constants", {}), overrides or {})
        self.materials = d.get("materials", {})
        if not isinstance(self.materials, dict):
            raise SystemError(f"{self.where}: 'materials' must be an object of name: material")
        for k in self.materials:
            self._check_name(k, "material")
        self.ambient = d.get("ambient", "Air")
        self.substrate = d["substrate"]
        self.layers = d["layers"]
        if not isinstance(self.layers, list):
            raise SystemError(f"{self.where}: 'layers' must be a list (light side first)")
        self._check_layers(self.layers, "layers")
        self.optics = self._optics(d.get("optics", {}))
        self.grids = d["sweep"] if isinstance(d["sweep"], list) else [d["sweep"]]
        self.rows = self._expand_sweep(self.grids)
        self.label = self._label(d.get("label", {}))
        self._cols: dict = {}
        self._nk_cache: dict = {}
        # resolve every material once up front so typos surface before the long calculation
        for i in self._sample_rows():
            self.stack(i)

    # -- validation helpers --------------------------------------------------------------------
    def _check_keys(self, obj: dict, allowed: set, where: str) -> None:
        bad = [k for k in obj if k not in allowed]
        if bad:
            raise SystemError(f"{self.where}: unknown key(s) {bad} in {where}; allowed: {sorted(allowed)}")

    def _check_name(self, name: str, what: str) -> None:
        if not _NAME.match(name) or name in _RESERVED:
            raise SystemError(f"{self.where}: '{name}' cannot be a {what} name (letters, digits and _, "
                              f"starting with a letter, and not a function name)")

    def _constants(self, raw: dict, overrides: dict) -> dict[str, float]:
        if not isinstance(raw, dict):
            raise SystemError(f"{self.where}: 'constants' must be an object of name: value")
        unknown = [k for k in overrides if k not in raw]
        if unknown:
            raise SystemError(f"{self.where}: cannot set {unknown}: not a constant of this system "
                              f"(constants: {', '.join(raw) or 'none'})")
        out: dict[str, float] = {}
        for k, v in raw.items():
            self._check_name(k, "constant")
            expr = overrides.get(k, v)
            try:
                val = float(np.asarray(evaluate(expr, out)))
            except (ExpressionError, TypeError) as e:
                raise SystemError(f"{self.where}: constant '{k}': {e}") from None
            out[k] = val
        return out

    def _optics(self, raw: dict) -> dict:
        self._check_keys(raw, _OPTICS, "'optics'")
        o = dict(illuminant="D65", na=0.0, na_weighting="uniform", angles=24)
        o.update({k: v for k, v in raw.items() if k != "notes"})
        o["na"] = float(evaluate(o["na"], self.constants))
        if not 0 <= o["na"] < 1:
            raise SystemError(f"{self.where}: optics.na must be in [0, 1) (objective numerical aperture in air)")
        if o["na_weighting"] not in ("uniform", "gaussian"):
            raise SystemError(f"{self.where}: optics.na_weighting must be 'uniform' or 'gaussian'")
        o["angles"] = int(o["angles"])
        return o

    def _check_layers(self, layers: list, where: str) -> None:
        for j, item in enumerate(layers):
            loc = f"{where}[{j}]" + (f" ('{item.get('name')}')" if isinstance(item, dict) and "name" in item else "")
            if not isinstance(item, dict):
                raise SystemError(f"{self.where}: {loc} must be an object with material and thickness")
            if "layers" in item:
                self._check_keys(item, _GROUP, loc)
                if not isinstance(item["layers"], list) or not item["layers"]:
                    raise SystemError(f"{self.where}: {loc}: 'layers' of a group must be a non-empty list")
                self._check_layers(item["layers"], loc + ".layers")
            else:
                self._check_keys(item, _LAYER, loc)
                for key in ("material", "thickness"):
                    if key not in item:
                        raise SystemError(f"{self.where}: {loc} needs '{key}'")

    # -- sweep -----------------------------------------------------------------------------------
    def _values(self, spec, pname: str) -> np.ndarray | None:
        """Swept values of one parameter, or None for a derived (expression) parameter."""
        if isinstance(spec, str):
            return None
        if isinstance(spec, bool):
            raise SystemError(f"{self.where}: sweep '{pname}': use numbers, not true/false")
        if isinstance(spec, (int, float)):
            return np.array([float(spec)])
        if isinstance(spec, list):
            return np.array([float(evaluate(v, self.constants)) for v in spec])
        if isinstance(spec, dict):
            ev = {k: float(evaluate(v, self.constants)) for k, v in spec.items() if k not in ("values", "notes")}
            if "values" in spec:
                return np.array([float(evaluate(v, self.constants)) for v in spec["values"]])
            if {"from", "to", "step"} <= set(spec):
                a, b, s = ev["from"], ev["to"], ev["step"]
                if s <= 0 or b < a:
                    raise SystemError(f"{self.where}: sweep '{pname}': need step > 0 and to >= from")
                n = int(np.floor((b - a) / s + 1e-9)) + 1
                return np.round(a + s * np.arange(n), 10)
            if {"from", "to", "num"} <= set(spec):
                return np.round(np.linspace(ev["from"], ev["to"], int(ev["num"])), 10)
        raise SystemError(f"{self.where}: sweep '{pname}': give a number, a list, {{\"from\", \"to\", \"step\"}}, "
                          f"{{\"from\", \"to\", \"num\"}}, {{\"values\": [...]}} or an expression string")

    def _expand_sweep(self, grids: list) -> dict[str, np.ndarray]:
        names = None
        parts = []
        for g, grid in enumerate(grids):
            if not isinstance(grid, dict) or not grid:
                raise SystemError(f"{self.where}: sweep[{g}] must be an object of parameter: values")
            grid = {k: v for k, v in grid.items() if k != "notes"}
            for k in grid:
                self._check_name(k, "sweep parameter")
                if k in self.constants:
                    raise SystemError(f"{self.where}: '{k}' is both a constant and a sweep parameter")
            if names is None:
                names = list(grid)
            elif set(grid) != set(names):
                raise SystemError(f"{self.where}: sweep[{g}] has parameters {sorted(grid)} but sweep[0] has "
                                  f"{sorted(names)}; every block must give every parameter")
            swept = {k: self._values(v, k) for k, v in grid.items()}
            fixed = [k for k, v in swept.items() if v is not None]
            if not fixed:
                raise SystemError(f"{self.where}: sweep[{g}] has no swept parameter (only expressions)")
            mesh = np.meshgrid(*[swept[k] for k in fixed], indexing="ij")  # first parameter varies slowest
            cols = {k: m.ravel() for k, m in zip(fixed, mesh)}
            for k, v in grid.items():  # derived parameters, in file order
                if swept[k] is None:
                    try:
                        cols[k] = np.broadcast_to(evaluate(v, {**self.constants, **cols}),
                                                  cols[fixed[0]].shape).astype(float)
                    except ExpressionError as e:
                        raise SystemError(f"{self.where}: sweep parameter '{k}': {e}") from None
            parts.append(cols)
        return {k: np.concatenate([p[k] for p in parts]) for k in names}

    def _sample_rows(self) -> list[int]:
        n = self.n_rows
        return sorted({0, n // 2, n - 1})

    @property
    def n_rows(self) -> int:
        return len(next(iter(self.rows.values())))

    @property
    def parameters(self) -> list[str]:
        """Sweep parameter names, in file order (a separately computed label column is not one)."""
        return self._grid_names()

    def _grid_names(self) -> list[str]:
        g = self.grids[0]
        return [k for k in g if k != "notes"]

    # -- label -----------------------------------------------------------------------------------
    def _label(self, raw: dict) -> dict:
        self._check_keys(raw, _LABEL, "'label'")
        params = self._grid_names()
        name = raw.get("name")
        if name is None:
            if len(params) != 1:
                raise SystemError(f"{self.where}: several sweep parameters ({params}); say which number the map "
                                  f"should report with label.name (and label.value if it is a combination)")
            name = params[0]
        self._check_name(name, "label")
        if name in params:
            if "value" in raw and str(raw["value"]).strip() != name:
                raise SystemError(f"{self.where}: label '{name}' is a sweep parameter; drop label.value or "
                                  f"give the label a new name")
        else:
            if "value" not in raw:
                raise SystemError(f"{self.where}: label '{name}' is not a sweep parameter, so give label.value, "
                                  f"an expression of {params}")
            try:
                val = evaluate(raw["value"], {**self.constants, **self.rows})
            except ExpressionError as e:
                raise SystemError(f"{self.where}: label.value: {e}") from None
            self.rows[name] = np.broadcast_to(np.asarray(val, float), (self.n_rows,)).copy()
        v = self.rows[name]
        classes = bool(raw.get("classes", False))
        span = float(np.nanmax(v) - np.nanmin(v)) or 1.0
        gap = raw.get("gap", 0.5 if classes else float(f"{span / 15:.3g}"))  # 3 significant figures
        return dict(name=name, unit=str(raw.get("unit", "")), title=str(raw.get("title", name)),
                    classes=classes, class_name=str(raw.get("class_name", "{}")),
                    zero_name=str(raw.get("zero_name", "")), gap=float(evaluate(gap, self.constants)))

    # -- expressions over rows -------------------------------------------------------------------
    def column(self, expr, what: str = "") -> np.ndarray:
        """Value of a number/expression for every candidate structure (cached)."""
        key = expr if isinstance(expr, str) else repr(expr)
        if key not in self._cols:
            try:
                val = evaluate(expr, {**self.constants, **self.rows})
            except ExpressionError as e:
                raise SystemError(f"{self.where}: {what}: {e}") from None
            except TypeError:
                raise SystemError(f"{self.where}: {what}: expected a number or expression, got {expr!r}") from None
            self._cols[key] = np.broadcast_to(np.asarray(val, dtype=float), (self.n_rows,))
        return self._cols[key]

    # -- materials -------------------------------------------------------------------------------
    def _file(self, name: str) -> Path | None:
        p = Path(name)
        cands = [p] if p.is_absolute() else ([self.path.parent / p] if self.path else []) + \
            [REPO / p, M.DATA_DIR / "materials" / p]
        return next((c for c in cands if c.is_file()), None)

    def _file_concrete(self, name: str, where: str):
        f = self._file(name)
        if f is None:
            raise SystemError(f"{self.where}: {where}: n,k file '{name}' not found (looked next to the system "
                              f"file, in the repository root and in data/materials/)")
        return ("file", str(f))

    def _concrete(self, spec, i: int, where: str, seen: tuple = ()):
        """The material spec for row i with every expression replaced by its number (hashable)."""
        if isinstance(spec, (int, float)) and not isinstance(spec, bool):
            return ("n", float(spec), 0.0)
        if isinstance(spec, str):
            if spec in self.materials:
                if spec in seen:
                    raise SystemError(f"{self.where}: material '{spec}' refers to itself")
                return self._concrete(self.materials[spec], i, f"material '{spec}'", seen + (spec,))
            if Path(spec).suffix.lower() in (".csv", ".txt", ".dat", ".nk", ".tsv"):
                return self._file_concrete(spec, where)
            if M._key(spec) not in M._library():
                raise SystemError(f"{self.where}: {where}: '{spec}' is not defined under 'materials' and not in "
                                  f"the n,k library ({', '.join(M.library_materials())})")
            return ("library", spec)
        if not isinstance(spec, dict):
            raise SystemError(f"{self.where}: {where}: a material is a name, a number or an object")
        spec = {k: v for k, v in spec.items() if k != "notes"}
        kind = next((k for k in _MATERIAL_KINDS if k in spec), None)
        if kind is None:
            raise SystemError(f"{self.where}: {where}: material object needs one of {sorted(_MATERIAL_KINDS)}")
        self._check_keys(spec, _MATERIAL_KINDS[kind], where)

        def num(key, default=None):
            if key not in spec:
                if default is None:
                    raise SystemError(f"{self.where}: {where}: missing '{key}'")
                return default
            return float(self.column(spec[key], f"{where}.{key}")[i])

        if kind == "library":
            return self._concrete(str(spec["library"]), i, where, seen)
        if kind == "file":
            return self._file_concrete(str(spec["file"]), where)
        if kind == "n":
            return ("n", num("n"), num("k", 0.0))
        if kind == "cauchy":
            c = spec["cauchy"]
            c = [c[k] for k in ("A", "B", "C") if k in c] if isinstance(c, dict) else c
            if not isinstance(c, list) or not 1 <= len(c) <= 3:
                raise SystemError(f"{self.where}: {where}: cauchy is [A, B, C] or {{\"A\":..,\"B\":..,\"C\":..}}")
            coeffs = tuple(float(self.column(v, f"{where}.cauchy")[i]) for v in c)
            return ("cauchy", coeffs, num("k", 0.0))
        if kind == "sellmeier":
            s = spec["sellmeier"]
            if not isinstance(s, dict) or set(s) != {"B", "C"}:
                raise SystemError(f"{self.where}: {where}: sellmeier is {{\"B\": [...], \"C\": [...]}} (C in um^2)")
            B = tuple(float(self.column(v, f"{where}.sellmeier.B")[i]) for v in s["B"])
            Cc = tuple(float(self.column(v, f"{where}.sellmeier.C")[i]) for v in s["C"])
            return ("sellmeier", B, Cc, num("k", 0.0))
        rule = str(spec["ema"]).lower()
        if rule not in M.EMA_RULES:
            raise SystemError(f"{self.where}: {where}: ema must be one of {sorted(M.EMA_RULES)}")
        for key in ("host", "inclusion", "fraction"):
            if key not in spec:
                raise SystemError(f"{self.where}: {where}: an ema material needs '{key}'")
        f = num("fraction")
        if not -1e-12 <= f <= 1 + 1e-12:
            raise SystemError(f"{self.where}: {where}: fraction = {f:g} for row {i} ({self._row_text(i)}); "
                              f"it must stay within [0, 1]")
        return ("ema", rule, self._concrete(spec["host"], i, where + ".host", seen),
                self._concrete(spec["inclusion"], i, where + ".inclusion", seen), min(max(f, 0.0), 1.0))

    def _nk(self, concrete, wl: np.ndarray) -> np.ndarray:
        if concrete in self._nk_cache:
            return self._nk_cache[concrete]
        kind = concrete[0]
        if kind == "n":
            out = np.full(wl.shape, complex(concrete[1], concrete[2]))
        elif kind in ("library", "file"):
            out = M.nk(concrete[1], wl)
        elif kind == "cauchy":
            out = M.cauchy(wl, *concrete[1]) + 1j * concrete[2]
        elif kind == "sellmeier":
            out = M.sellmeier(wl, list(concrete[1]), list(concrete[2])).real + 1j * concrete[3]
        else:
            _, rule, host, incl, f = concrete
            out = M.EMA_RULES[rule](self._nk(host, wl), self._nk(incl, wl), f)
        self._nk_cache[concrete] = out
        return out

    # -- stacks ----------------------------------------------------------------------------------
    def _row_text(self, i: int) -> str:
        return ", ".join(f"{k}={v[i]:g}" for k, v in self.rows.items())

    def _expand(self, layers: list, i: int, where: str) -> list[tuple[dict, float, str]]:
        out = []
        for j, item in enumerate(layers):
            loc = f"{where}[{j}]"
            name = item.get("name", loc)
            loc = f"{loc} ('{name}')" if "name" in item else loc
            rep = float(self.column(item.get("repeat", 1), f"{loc} repeat")[i])
            if rep < -1e-9 or abs(rep - round(rep)) > 1e-6:
                raise SystemError(f"{self.where}: {loc}: repeat = {rep:g} for {self._row_text(i)}; "
                                  f"it must be a whole number >= 0")
            rep = int(round(rep))
            if "layers" in item:
                for _ in range(rep):
                    out += self._expand(item["layers"], i, f"{where}[{j}].layers")
                continue
            d = float(self.column(item["thickness"], f"{loc} thickness")[i])
            if not np.isfinite(d) or d < 0:
                raise SystemError(f"{self.where}: {loc}: thickness = {d:g} nm for {self._row_text(i)}")
            if d > 0:  # a zero-thickness layer changes nothing (I_ab I_bc = I_ac), so it is left out
                out += [(item, d, name)] * rep
        return out

    def stack(self, i: int, wl: np.ndarray | None = None):
        """(complex indices, thicknesses nm, names) of candidate ``i``, ambient first, substrate last."""
        from .colorimetry import WAVELENGTHS
        wl = WAVELENGTHS if wl is None else wl
        layers = self._expand(self.layers, i, "layers")
        N = [self._nk(self._concrete(self.ambient, i, "ambient"), wl)]
        N += [self._nk(self._concrete(item["material"], i, f"layer '{name}'"), wl) for item, _, name in layers]
        N.append(self._nk(self._concrete(self.substrate, i, "substrate"), wl))
        d = [0.0] + [t for _, t, _ in layers] + [0.0]
        names = ["ambient"] + [name for _, _, name in layers] + ["substrate"]
        return N, d, names

    def spectra(self, na: float | None = None) -> np.ndarray:
        """Reflectance of every candidate on 360-830 nm (rows x 471)."""
        from .colorimetry import WAVELENGTHS
        na = self.optics["na"] if na is None else na
        out = np.empty((self.n_rows, len(WAVELENGTHS)))
        for i in range(self.n_rows):
            N, d, _ = self.stack(i)
            out[i] = reflectance_na(N, d, WAVELENGTHS, na, self.optics["angles"], self.optics["na_weighting"])
        return out

    def at(self, assignments: list[dict[str, float]]) -> "System":
        """A copy whose candidates are the given parameter sets, e.g. [{'thickness_nm': 100}]."""
        grid0 = {k: v for k, v in self.grids[0].items() if k != "notes"}
        derived = {k: v for k, v in grid0.items() if isinstance(v, str)}
        needed = [k for k in grid0 if k not in derived]
        grids = []
        for a in assignments:
            missing = [k for k in needed if k not in a]
            extra = [k for k in a if k not in grid0]
            if missing or extra:
                raise SystemError(f"{self.where}: give exactly the parameters {needed} (missing {missing}, "
                                  f"unknown {extra})")
            grids.append({k: (float(a[k]) if k in a else derived[k]) for k in grid0})
        d = copy.deepcopy(self.definition)
        d["sweep"] = grids
        d["constants"] = {**d.get("constants", {}), **self.constants}
        if "gap" not in d.get("label", {}):
            d.setdefault("label", {})["gap"] = self.label["gap"]
        return System(d, self.path)

    # -- description -----------------------------------------------------------------------------
    def _material_text(self, spec) -> str:
        if isinstance(spec, (int, float)):
            return f"n={spec:g}"
        if isinstance(spec, str):
            return spec
        if "library" in spec:
            return str(spec["library"])
        if "file" in spec:
            return Path(str(spec["file"])).name
        if "n" in spec:
            return f"n={spec['n']}" + (f"+{spec['k']}i" if spec.get("k") not in (None, 0) else "")
        if "cauchy" in spec:
            return "Cauchy"
        if "sellmeier" in spec:
            return "Sellmeier"
        return (f"{spec['ema']}({self._material_text(spec['inclusion'])} in {self._material_text(spec['host'])}, "
                f"f={self._value_text(spec['fraction'])})")

    def _value_text(self, v) -> str:
        if isinstance(v, str) and not (names_in(v) - set(self.constants)):
            v = float(evaluate(v, self.constants))
        return f"{v:.4g}" if isinstance(v, (int, float)) else v

    def describe(self) -> str:
        """One-line stack, e.g. 'Air | flake: aMoO3 thickness_nm nm | oxide: SiO2-Franta 100 nm | Si-Franta'."""
        def walk(layers):
            parts = []
            for item in layers:
                rep = item.get("repeat", 1)
                rep_txt = "" if rep in (1, "1") else f" (repeat {self._value_text(rep)})"
                if "layers" in item:
                    parts.append(f"[{' | '.join(walk(item['layers']))}]{rep_txt}")
                else:
                    name = item.get("name", "")
                    thick = self._value_text(item["thickness"])
                    thick = f"({thick}) nm" if " " in thick else f"{thick} nm"
                    parts.append(f"{name + ': ' if name else ''}{self._material_text(item['material'])}, {thick}"
                                 f"{rep_txt}")
            return parts
        return " | ".join([self._material_text(self.ambient)] + walk(self.layers) +
                          [self._material_text(self.substrate)])

    def relative_path(self) -> str:
        if not self.path:
            return ""
        try:
            return self.path.resolve().relative_to(REPO).as_posix()
        except ValueError:
            return str(self.path)
