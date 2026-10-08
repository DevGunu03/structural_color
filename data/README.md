# Data shared by the Python and MATLAB code

| File | Content |
| --- | --- |
| `nk_library.csv` | Optical constants n and k of 31 materials, 300–800 nm in 1 nm steps: a CSV export of `Index_of_Refraction_library.xls` used by the original transfer-matrix code (`legacy/TransferMatrix/`) |
| `cie1931_2deg.csv` | CIE 1931 2° colour-matching functions x̄, ȳ, z̄, 360–830 nm in 1 nm steps |
| `cie_illuminants.csv` | relative spectral power of CIE illuminants D65, A and D50, 360–830 nm |
| `materials/` | materials added by users, one file each, usable by name like the built-in ones ([how to add one](materials/README.md)) |

## The n,k library

The library is `nk_library.csv` plus the files in `materials/`. In `nk_library.csv`, each material is a pair of columns `<name>_n`, `<name>_k`. In system files and on the command line, case, `-`, `_` and spaces in names are ignored, so `SiO2-Franta` and `sio2_franta` both work. List the names with `python python/run_tcd.py materials` or `tcd_materials` in MATLAB. Print one material's values with `python python/run_tcd.py materials --show aMoO3`.

The entries used by the bundled systems come from:

- `Si_Franta`, `SiO2_Franta`, `TiO2_Franta`: Franta et al.
- `PS_beads`: Cauchy fit, Naglič et al. 2020.
- `aMoO3`: α-MoO₃, Lajaunie et al.

The other entries (metals, organic semiconductors, other TiO₂ films, ...) were already in the McGehee group's library or were added during earlier projects. Check their source before relying on them.

Values between tabulated wavelengths are interpolated linearly. Outside the table they are **extrapolated** linearly, as in the original code. The library stops at 800 nm, so 800–830 nm, where the eye is nearly blind, is always extrapolated.

### Adding a material

Use the `add-material` command (Python) or `tcd_add_material` (MATLAB). The command examines the data and saves them in `data/materials/<name>.csv`. Then all system files can use the material by name, for example `"material": "MoS2"`. For the procedure, refer to [materials/README.md](materials/README.md).

Do not add columns to `nk_library.csv` by hand. The built-in library stops at 800 nm, and a file in `data/materials/` keeps the source of the data.

Optical constants and the oxide thickness are the largest causes of errors in the calculated colours. Always record the source of the data.

## Colorimetry tables

Both tables were exported from [colour-science](https://www.colour-science.org/) 0.4.6 (`MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]`, `SDS_ILLUMINANTS`) so that MATLAB, which has no built-in CIE data, uses identical numbers. A measured lamp spectrum can be used instead of these illuminants: give a two-column CSV (wavelength nm, power) as the illuminant.
