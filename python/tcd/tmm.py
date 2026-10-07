"""Transfer-matrix reflectance of a planar multilayer.

Same formalism and sign conventions as the McGehee-group TransferMatrix.m
(Pettersson et al., J. Appl. Phys. 86, 487 (1999)) used in TransferMatrix_Updated/:
interface matrix I = [[1, r], [r, 1]] / t, layer matrix L = diag(exp(-i xi d), exp(i xi d)),
S = I_01 L_1 I_12 ... L_{N-1} I_{N-1,N}, r = S21 / S11. The first medium (air) is
the incidence side and the last medium (Si) is semi-infinite.

At normal incidence this reproduces the MATLAB output to machine precision
(see tests/validate_against_matlab.py). ``reflectance_na`` adds the oblique s/p
calculation integrated over the objective's numerical aperture, which the
MATLAB pipeline did not include.
"""
from __future__ import annotations

import numpy as np


def _interface(n1, n2, c1=None, c2=None, pol="s"):
    if c1 is None:  # normal incidence
        r = (n1 - n2) / (n1 + n2)
        t = 2 * n1 / (n1 + n2)
    elif pol == "s":
        a, b = n1 * c1, n2 * c2
        r = (a - b) / (a + b)
        t = 2 * a / (a + b)
    else:  # p polarisation
        a, b = n2 * c1, n1 * c2
        r = (a - b) / (a + b)
        t = 2 * n1 * c1 / (a + b)
    m = np.empty(r.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = 1 / t
    m[..., 0, 1] = r / t
    m[..., 1, 0] = r / t
    m[..., 1, 1] = 1 / t
    return m


def _layer(phase):
    m = np.zeros(phase.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = np.exp(-1j * phase)
    m[..., 1, 1] = np.exp(1j * phase)
    return m


def reflectance(n_stack: list[np.ndarray], thicknesses: list[float], wavelengths: np.ndarray) -> np.ndarray:
    """Normal-incidence power reflectance R(lambda).

    n_stack     : complex index arrays, one per medium, incidence medium first, substrate last.
    thicknesses : nm, same length as n_stack; first and last entries are ignored.
    """
    wl = np.asarray(wavelengths, dtype=float)
    S = _interface(n_stack[0], n_stack[1])
    for j in range(1, len(n_stack) - 1):
        S = S @ _layer(2 * np.pi * n_stack[j] * thicknesses[j] / wl) @ _interface(n_stack[j], n_stack[j + 1])
    return np.abs(S[..., 1, 0] / S[..., 0, 0]) ** 2


def reflectance_oblique(n_stack, thicknesses, wavelengths, theta0, pol):
    """Reflectance at incidence angle theta0 (rad, in the first medium) for 's' or 'p'."""
    wl = np.asarray(wavelengths, dtype=float)
    kx = n_stack[0] * np.sin(theta0)  # conserved tangential component
    cos = []
    for n in n_stack:
        c = np.sqrt(1 - (kx / n) ** 2 + 0j)
        c = np.where(np.imag(n * c) < 0, -c, c)  # forward-decaying branch
        cos.append(c)
    S = _interface(n_stack[0], n_stack[1], cos[0], cos[1], pol)
    for j in range(1, len(n_stack) - 1):
        phase = 2 * np.pi * n_stack[j] * cos[j] * thicknesses[j] / wl
        S = S @ _layer(phase) @ _interface(n_stack[j], n_stack[j + 1], cos[j], cos[j + 1], pol)
    return np.abs(S[..., 1, 0] / S[..., 0, 0]) ** 2


def reflectance_na(n_stack, thicknesses, wavelengths, na: float, n_angles: int = 24,
                   weighting: str = "uniform") -> np.ndarray:
    """Unpolarised reflectance averaged over the illumination cone of an objective.

    R = int (Rs + Rp)/2 w(theta) sin(theta) cos(theta) dtheta / int w(theta) sin(theta) cos(theta) dtheta,
    theta in [0, asin(NA)]. ``weighting`` is 'uniform' (filled back aperture) or
    'gaussian' (w = exp(-theta^2 / theta_b^2), theta_b = asin(NA), as in Jung et al. 2007).
    """
    if na <= 0:
        return reflectance(n_stack, thicknesses, wavelengths)
    tmax = np.arcsin(na)
    x, wq = np.polynomial.legendre.leggauss(n_angles)
    theta = 0.5 * tmax * (x + 1)
    wq = 0.5 * tmax * wq
    w = np.exp(-(theta / tmax) ** 2) if weighting == "gaussian" else np.ones_like(theta)
    w = w * wq * np.sin(theta) * np.cos(theta)
    R = np.zeros(np.shape(wavelengths))
    for th, wi in zip(theta, w):
        R += wi * 0.5 * (reflectance_oblique(n_stack, thicknesses, wavelengths, th, "s")
                         + reflectance_oblique(n_stack, thicknesses, wavelengths, th, "p"))
    return R / w.sum()
