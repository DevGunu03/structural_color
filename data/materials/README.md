# Added materials

This folder contains optical constants (n and k) that users added to the library. Each file is one material. The file name is the material name.

All system files can use these materials by name, in Python and in MATLAB:

```jsonc
"layers": [{"name": "flake", "material": "MoS2", "thickness": "t"}]
```

## Add a material

1. Get a file with the wavelength, n and k of your material. Refer to [Accepted files](#accepted-files).
2. Type this command in the repository folder. Change the name, the file and the source.

   ```bash
   python python/run_tcd.py add-material MoS2 MoS2_downloaded.csv --source "Author et al. 2020, refractiveindex.info"
   ```

   In MATLAB, in the `matlab` folder:

   ```matlab
   tcd_add_material('MoS2', 'MoS2_downloaded.csv', 'Source', 'Author et al. 2020, refractiveindex.info')
   ```

3. Read the result. The software shows n and k at 450 nm, 550 nm and 650 nm.
4. Read the warnings. A warning tells you if the data do not cover 360–830 nm.
5. Type this command to see the new material in the library:

   ```bash
   python python/run_tcd.py materials
   ```

You can also give these options:

| Option (Python) | Option (MATLAB) | Use |
| --- | --- | --- |
| `--reference` | `'Reference'` | A DOI or a URL for the data. |
| `--notes` | `'Notes'` | Other information, for example the crystal axis, or "film" or "bulk". |
| `--units nm` or `--units um` | `'Units'` | The wavelength unit, if the automatic choice is wrong. |
| `--replace` | `'Replace', true` | Replace an added material that has the same name. |

The software does not add the material if:

- The name is already in the built-in library (`data/nk_library.csv`).
- The name is already in this folder, and you did not give `--replace`.
- The name has characters other than letters, digits, `-` and `_`.
- You did not give a source.
- n is zero or less, or k is less than zero.

## Accepted files

The software reads these types of file:

- **A table** with the columns wavelength, n and k. The columns can have commas, tabs, spaces or semicolons between them. If there is no k column, k is zero.
- **A CSV file from [refractiveindex.info](https://refractiveindex.info)**. This file has a table of n, then a table of k. The software joins the two tables.

Header lines and comment lines are not a problem. The software ignores them.

The wavelength can be in nm or in µm. If all wavelengths are less than 50, the software uses µm.

> **Note:** The colours use the wavelengths from 360 nm to 830 nm. Outside the range of your data, the software extrapolates n and k in a straight line. Near an absorption edge, this extrapolation can be very wrong. Use data that cover at least 400–700 nm.

## The file format

The software writes each material in this format:

```text
# name: MyMaterial
# source: Author et al. 2020, refractiveindex.info
# reference: https://doi.org/...
# notes: example values only
# added: 2026-10-08
# imported from: MyMaterial_downloaded.csv
wavelength_nm,n,k
360,2.41,0.05
370,2.39,0.04
...
```

You can also write a file in this format by hand. Put it in this folder. The software finds it automatically.

## Share a material

Other users can use your material if you add it to the repository.

1. Add the material with the `add-material` command.
2. Make sure that the source is correct and complete.
3. Commit the new file in `data/materials/`.
4. Open a pull request on GitHub. In the description, tell the reviewers where the data come from.

> **Note:** Do not share data that you cannot share legally. Some publications do not allow it. Data from refractiveindex.info show their license on the page of each material.
