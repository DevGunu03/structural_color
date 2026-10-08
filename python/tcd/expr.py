"""Arithmetic expressions used in system files, e.g. ``"thickness": "bead_nm * sqrt(2/3)"``.

The expression is parsed by a small recursive-descent parser and evaluated on numbers or
numpy arrays (one element per candidate structure). Nothing is passed to Python's ``eval``,
so a system file cannot run code. The MATLAB twin is ``matlab/private/expr_eval.m``; both
accept exactly the same language:

    numbers      2, 0.5, .5, 1e-3
    names        constants and sweep parameters of the system file, and ``pi``
    operators    + - * / ^ (power, ** also accepted), unary -, parentheses
    comparisons  < <= > >= == !=   (give 1 for true, 0 for false)
    functions    sqrt exp log log10 sin cos tan asin acos atan abs floor ceil round
                 min(a, b) max(a, b)

Precedence, low to high: comparison, + -, * /, unary -, ^ (right-associative), so
``-2^2 = -4`` and ``2^-1 = 0.5`` as in MATLAB and Python. ``round`` rounds halves up
(``round(2.5) = 3``, ``round(-2.5) = -2``) in both languages.
"""
from __future__ import annotations

import re

import numpy as np

_TOKEN = re.compile(r"\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)"  # number
                    r"|([A-Za-z_]\w*)"                                           # name
                    r"|(\*\*|<=|>=|==|!=|[-+*/^(),<>]))")                        # operator
_FUNCS = {
    "sqrt": np.sqrt, "exp": np.exp, "log": np.log, "log10": np.log10,
    "sin": np.sin, "cos": np.cos, "tan": np.tan, "asin": np.arcsin, "acos": np.arccos, "atan": np.arctan,
    "abs": np.abs, "floor": np.floor, "ceil": np.ceil, "round": lambda x: np.floor(np.asarray(x) + 0.5),
    "min": np.minimum, "max": np.maximum,
}
_NARGS = {"min": 2, "max": 2}
_CMP = {"<": np.less, "<=": np.less_equal, ">": np.greater, ">=": np.greater_equal,
        "==": np.equal, "!=": np.not_equal}


class ExpressionError(ValueError):
    pass


def _tokens(text: str) -> list[tuple[str, str]]:
    out, pos, text = [], 0, text.rstrip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ExpressionError(f"cannot read '{text[pos:].strip()}' in expression '{text}'")
        num, name, op = m.groups()
        out.append(("num", num) if num else ("name", name) if name else ("op", "^" if op == "**" else op))
        pos = m.end()
    return out


def names_in(text: str) -> set[str]:
    """Variable names an expression refers to (function names excluded)."""
    toks = _tokens(str(text))
    return {v for i, (k, v) in enumerate(toks) if k == "name" and v != "pi"
            and not (v in _FUNCS and i + 1 < len(toks) and toks[i + 1] == ("op", "("))}


def evaluate(text, variables: dict | None = None):
    """Value of ``text`` (a number or an expression string) with ``variables`` substituted.

    Variables may be numbers or equal-length numpy arrays; the result broadcasts like numpy.
    """
    if isinstance(text, (int, float, np.number, np.ndarray)) and not isinstance(text, bool):
        return text
    if not isinstance(text, str):
        raise ExpressionError(f"expected a number or an expression string, got {text!r}")
    variables = variables or {}
    toks = _tokens(text)
    pos = 0

    def peek():
        return toks[pos] if pos < len(toks) else ("end", "")

    def take(expected=None):
        nonlocal pos
        tok = peek()
        if expected is not None and tok != ("op", expected):
            got = tok[1] or "end of expression"
            raise ExpressionError(f"expected '{expected}' but found '{got}' in '{text}'")
        pos += 1
        return tok

    def comparison():
        left = additive()
        if peek()[0] == "op" and peek()[1] in _CMP:
            op = take()[1]
            left = _CMP[op](left, additive()).astype(float)
        return left

    def additive():
        val = term()
        while peek() in (("op", "+"), ("op", "-")):
            op = take()[1]
            rhs = term()
            val = val + rhs if op == "+" else val - rhs
        return val

    def term():
        val = unary()
        while peek() in (("op", "*"), ("op", "/")):
            op = take()[1]
            rhs = unary()
            with np.errstate(divide="ignore", invalid="ignore"):
                val = val * rhs if op == "*" else np.true_divide(val, rhs)
        return val

    def unary():
        if peek() in (("op", "-"), ("op", "+")):
            sign = -1.0 if take()[1] == "-" else 1.0
            return sign * unary()
        return power()

    def power():
        base = atom()
        if peek() == ("op", "^"):
            take()
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.power(np.asarray(base, dtype=float), unary())
        return base

    def atom():
        kind, val = take()
        if kind == "num":
            return float(val)
        if kind == "name":
            if peek() == ("op", "(") and val in _FUNCS:
                take("(")
                args = [comparison()]
                while peek() == ("op", ","):
                    take()
                    args.append(comparison())
                take(")")
                need = _NARGS.get(val, 1)
                if len(args) != need:
                    raise ExpressionError(f"{val}() takes {need} argument(s), got {len(args)} in '{text}'")
                with np.errstate(divide="ignore", invalid="ignore"):
                    return _FUNCS[val](*args)
            if val == "pi" and "pi" not in variables:
                return np.pi
            if val not in variables:
                known = ", ".join(sorted(variables)) or "none"
                raise ExpressionError(f"unknown name '{val}' in '{text}' (known names: {known})")
            return variables[val]
        if (kind, val) == ("op", "("):
            inner = comparison()
            take(")")
            return inner
        raise ExpressionError(f"unexpected '{val or 'end of expression'}' in '{text}'")

    if not toks:
        raise ExpressionError("empty expression")
    result = comparison()
    if pos != len(toks):
        raise ExpressionError(f"unexpected '{toks[pos][1]}' in '{text}'")
    return result
