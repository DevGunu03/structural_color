# Your own optical constants

Put (wavelength, n, k) files here and use them in a system file by path:

```jsonc
"layers": [{"name": "flake", "material": "data/materials/MoS2.csv", "thickness": "t"}]
```

A file name is resolved next to the system file first, then from the repository root, then in this folder. `"material": "MoS2.csv"` therefore also finds `data/materials/MoS2.csv`.

## Format

```text
# MoS2, bulk, from <reference / DOI / your ellipsometry fit>
wavelength_nm, n, k
350, 3.10, 2.21
360, 3.25, 2.05
...
```

- **Columns.** Wavelength, n, and optionally k (k = 0 if absent). k > 0 means absorption. Extra columns are ignored.
- **Separators.** Commas, tabs, spaces or semicolons. Lines that do not start with numbers (headers, comments) are skipped.
- **Units.** Wavelength in nm, or in µm if every value is below 50 (as exported by refractiveindex.info).
- **Range.** Cover 360–830 nm if you can. Outside the tabulated range the values are extrapolated linearly, which can go badly wrong near absorption edges.
- **File name.** It must end in `.csv`, `.txt`, `.dat`, `.nk` or `.tsv`.

Check the file before building a reference:

```bash
python python/run_tcd.py materials --show data/materials/MoS2.csv
```
```matlab
n = tcd_materials('../data/materials/MoS2.csv', 400:50:800)
```
