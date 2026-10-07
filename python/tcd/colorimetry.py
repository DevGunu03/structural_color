"""Reflectance spectra -> CIE XYZ -> CAM02-UCS / CIELAB / sRGB, and back from camera sRGB.

Conventions (kept identical to the existing pipeline so results are comparable):
  * CIE 1931 2-degree observer, 360-830 nm at 1 nm.
  * XYZ scaled so a perfect reflector has Y = 100 under the chosen illuminant.
  * CAM02-UCS viewing conditions as in sRGB_to_CAM02UCS.m (Cobeldick CIECAM02 toolbox):
    D65 white, L_A = 64 / (pi * 5) cd/m^2, Y_b = 20, average surround.

Fixes relative to the notebooks:
  * CIELAB uses the illuminant's xy chromaticity as the white point. The notebooks
    passed sd_to_XYZ(illuminant, cmfs, illuminant) (Y ~ 10^4) as the white, which
    produced L* ~ 1-3 for every sample.
  * CAM02-UCS is computed straight from XYZ, without clipping to the sRGB gamut first.
"""
from __future__ import annotations

from pathlib import Path

import colour
import numpy as np

SHAPE = colour.SpectralShape(360, 830, 1)
WAVELENGTHS = np.arange(360, 831, 1, dtype=float)
CMFS = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(SHAPE)
XY_D65 = colour.CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["D65"]
XYZ_D65 = colour.xy_to_XYZ(XY_D65) * 100

CAM02_VIEWING = dict(XYZ_w=XYZ_D65, L_A=64 / np.pi / 5, Y_b=20.0,
                     surround=colour.VIEWING_CONDITIONS_CIECAM02["Average"])


def illuminant_spd(illuminant: str | Path = "D65") -> np.ndarray:
    """Relative spectral power on WAVELENGTHS. Accepts a CIE name ('D65', 'A', 'D50', ...)
    or a two-column CSV (wavelength nm, power), e.g. a measured microscope-lamp spectrum."""
    name = str(illuminant)
    if name in colour.SDS_ILLUMINANTS:
        return colour.SDS_ILLUMINANTS[name].copy().align(SHAPE).values
    data = np.genfromtxt(name, delimiter=",")
    data = data[~np.isnan(data).any(axis=1)]
    return np.interp(WAVELENGTHS, data[:, 0], data[:, 1], left=0.0, right=0.0)


def reflectance_to_XYZ(R: np.ndarray, illuminant: str | Path = "D65", adapt_to_d65: bool = True) -> np.ndarray:
    """Tristimulus values of reflectance spectra R[..., 471] sampled on WAVELENGTHS.

    With a non-D65 illuminant and ``adapt_to_d65`` the result is chromatically adapted
    (CAT02) to D65, which is what a white-balanced camera does to the scene.
    """
    R = np.asarray(R, dtype=float)
    S = illuminant_spd(illuminant)
    xbar, ybar, zbar = CMFS.values.T
    k = 100.0 / np.sum(S * ybar)
    XYZ = k * np.stack([np.sum(R * S * xbar, axis=-1),
                        np.sum(R * S * ybar, axis=-1),
                        np.sum(R * S * zbar, axis=-1)], axis=-1)
    if adapt_to_d65 and str(illuminant) != "D65":
        white = k * np.array([np.sum(S * xbar), np.sum(S * ybar), np.sum(S * zbar)])
        XYZ = colour.chromatic_adaptation(XYZ, white, XYZ_D65, method="Von Kries", transform="CAT02")
    return XYZ


def XYZ_to_cam02ucs(XYZ: np.ndarray) -> np.ndarray:
    """J', a', b' (CAM02-UCS) for XYZ on the 0-100 scale, D65-referred.

    Goes through XYZ_to_CIECAM02 explicitly: colour.XYZ_to_CAM02UCS expects XYZ on a
    0-1 scale, and feeding it 0-100 values silently gives J' ~ 25 % too high.
    """
    spec = colour.XYZ_to_CIECAM02(np.asarray(XYZ, dtype=float), CAM02_VIEWING["XYZ_w"],
                                  CAM02_VIEWING["L_A"], CAM02_VIEWING["Y_b"], CAM02_VIEWING["surround"])
    return colour.JMh_CIECAM02_to_CAM02UCS(np.stack([spec.J, spec.M, spec.h], axis=-1))


def XYZ_to_lab(XYZ: np.ndarray) -> np.ndarray:
    return colour.XYZ_to_Lab(np.asarray(XYZ, dtype=float) / 100.0, XY_D65)


def XYZ_to_srgb(XYZ: np.ndarray, clip: bool = True) -> np.ndarray:
    rgb = colour.XYZ_to_sRGB(np.asarray(XYZ, dtype=float) / 100.0)
    return np.clip(rgb, 0, 1) if clip else rgb


# --- camera side -------------------------------------------------------------------
SRGB = colour.RGB_COLOURSPACES["sRGB"]


def srgb_to_linear(rgb: np.ndarray) -> np.ndarray:
    return colour.models.eotf_sRGB(np.asarray(rgb, dtype=float))


def linear_to_XYZ(rgb_linear: np.ndarray) -> np.ndarray:
    """Linear sRGB (D65) -> XYZ on the 0-100 scale."""
    return 100.0 * np.einsum("ij,...j->...i", SRGB.matrix_RGB_to_XYZ, rgb_linear)


def XYZ_to_linear(XYZ: np.ndarray) -> np.ndarray:
    return np.einsum("ij,...j->...i", SRGB.matrix_XYZ_to_RGB, np.asarray(XYZ, dtype=float) / 100.0)


def delta_e(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Euclidean colour difference (Delta E in CAM02-UCS, or Delta E*ab in CIELAB)."""
    return np.linalg.norm(np.asarray(a) - np.asarray(b), axis=-1)
