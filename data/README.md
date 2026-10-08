# Data shared by the Python and MATLAB code

| File | Content |
| --- | --- |
| `nk_library.csv` | Optical constants n and k of 31 materials, 300–800 nm in 1 nm steps: a CSV export of `Index_of_Refraction_library.xls` used by the original transfer-matrix code (`legacy/TransferMatrix/`) |
| `cie1931_2deg.csv` | CIE 1931 2° colour-matching functions x̄, ȳ, z̄, 360–830 nm in 1 nm steps |
| `cie_illuminants.csv` | relative spectral power of CIE illuminants D65, A and D50, 360–830 nm |
| `materials/` | your own (wavelength, n, k) files, usable by file name in system files ([format](materials/README.md)) |

## The n,k library

Each material is a pair of columns `<name>_n`, `<name>_k`. In system files and on the command line, case, `-`, `_` and spaces in names are ignored, so `SiO2-Franta` and `sio2_franta` both work. List the names with `python python/run_tcd.py materials` or `tcd_materials` in MATLAB. Print one material's values with `python python/run_tcd.py materials --show aMoO3`.

The entries used by the bundled systems come from:

- `Si_Franta`, `SiO2_Franta`, `TiO2_Franta`: Franta et al.
- `PS_beads`: Cauchy fit, Naglič et al. 2020.
- `aMoO3`: α-MoO₃, Lajaunie et al.

The other entries (metals, organic semiconductors, other TiO₂ films, ...) were already in the McGehee group's library or were added during earlier projects. Check their source before relying on them.

Values between tabulated wavelengths are interpolated linearly. Outside the table they are **extrapolated** linearly, as in the original code. The library stops at 800 nm, so 800–830 nm, where the eye is nearly blind, is always extrapolated.

### Adding a material

There are two ways:

1. **A file** in `data/materials/` (recommended). Write `"material": "data/materials/MoS2.csv"` in the system file. Nothing else changes.
2. **A column pair in `nk_library.csv`** (`MoS2_n`, `MoS2_k`), interpolated onto the library's 300–800 nm grid. The name then works everywhere, like the built-in materials. `legacy/TransferMatrix/addNK_to_TransferMatrix.m` did this for the old Excel library.

Note where the data came from in your system file's comments. Optical constants are the largest single source of error in the simulated colours.

## Colorimetry tables

Both tables were exported from [colour-science](https://www.colour-science.org/) 0.4.6 (`MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]`, `SDS_ILLUMINANTS`) so that MATLAB, which has no built-in CIE data, uses identical numbers. A measured lamp spectrum can be used instead of these illuminants: give a two-column CSV (wavelength nm, power) as the illuminant.
