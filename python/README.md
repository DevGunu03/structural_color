# Python: `run_tcd.py` and the `tcd` package

Python 3.10 or newer. Install the dependencies once, then run everything from the repository root:

```bash
pip install -r python/requirements.txt
python python/tests/test_pipeline.py      # prints "all checks passed"
```

## Commands

`python python/run_tcd.py <command> --help` lists every option.

### `materials`: what is available

```bash
python python/run_tcd.py materials                 # n,k library names and bundled systems
python python/run_tcd.py materials --show aMoO3    # n and k at 400-800 nm (a library name or a file)
```

### `spectrum`: check single structures of a system file

```bash
python python/run_tcd.py spectrum --system ps --at layers=1,packing=1 --at layers=2,packing=0.5
python python/run_tcd.py spectrum --system systems/my_sample.jsonc --at t=100 --csv
```

This prints each stack as built, with its XYZ, J′a′b′, sRGB and ΔE from the first structure, and writes a plot of R(λ) with colour swatches (`results/<name>_spectra.png`, plus a `.csv` with `--csv`). Give every sweep parameter in `--at`; derived parameters are computed. Without `--at`, it shows the first, middle and last candidates.

### `build-ref`: simulate a system's colour reference

```bash
python python/run_tcd.py build-ref --system moo3                              # -> refs/moo3_D65.csv
python python/run_tcd.py build-ref --system ps --set bead_nm=500 --set oxide_nm=285
python python/run_tcd.py build-ref --system systems/my_sample.jsonc --na 0.5 --spectra
```

| Option | Meaning |
| --- | --- |
| `--system` | a system file, or the name of one in `systems/` |
| `--set NAME=VALUE` | change a constant of the system file (repeatable) |
| `--illuminant` | override `optics.illuminant`: `D65`, `A`, `D50` or a lamp-spectrum CSV |
| `--na` | override `optics.na` (objective NA; 0 = normal incidence) |
| `--out` | reference CSV; default `refs/<name>[_<changes>][_NA<na>]_<illuminant>.csv` |
| `--spectra` | also save every reflectance spectrum (`<out>_spectra.csv`) |

Each run writes:

- `<out>.csv`: one row per candidate, with the sweep parameters, the label, XYZ, J′a′b′, L\*a\*b\*, sRGB and ΔE from row 0;
- `<out>.json`: the label description, the stack, the constants and the full system definition;
- `<out>.png`: the reference sheet (see [systems/README.md](../systems/README.md#4-build-the-reference-and-read-its-sheet)).

### `map`: thickness / layer map of a micrograph

```bash
python python/run_tcd.py map --ref refs/ps_D65.csv --interactive                    # choose the image, click the substrate
python python/run_tcd.py map --ref refs/moo3_D65.csv --image flakes.png --substrate auto --max-value 350
python python/run_tcd.py map --ref refs/ps_D65.csv --image beads.png --substrate-colour 0.157,0.263,0.459 \
                             --anchor 1L@120,40,30,30 --out results/beads
```

| Option | Meaning |
| --- | --- |
| `--ref` | reference CSV from `build-ref` (any system) |
| `--image` | the micrograph; without it, a file dialog opens |
| `--substrate` | bare-substrate region `x,y,w,h` (pixels, 0-based top-left), `auto` (dominant colour) or `colour:r,g,b` |
| `--substrate-colour r,g,b` | find the substrate region from its colour (0–1) |
| `--anchor LABEL@REGION` | an extra region of known structure: `1L@x,y,w,h`, `250nm@colour:r,g,b`, `layers=2,packing=1@x,y,w,h` (repeatable) |
| `--interactive` | click any region not given on the command line |
| `--gamma` | camera decoding: a number g (intensity = value^g; default 1, a linear camera) or `srgb` |
| `--max-value` | only consider candidates whose label is ≤ this, e.g. a thickness bound from AFM (`--t-max` and `--max-layers` are older names) |
| `--gap` | ambiguity gap in label units (default from the reference) |
| `--max-residual` | ΔE above which a pixel stays unassigned (default 15) |
| `--correction` | `auto` (gains for one region, 3×3 matrix for more), `diagonal` or `matrix` |
| `--sigma`, `--max-side` | Gaussian smoothing (px, default 1.2) and downsampling of the longer side (default 1600; 0 = off) |

Each run writes:

- `<out>.png`: micrograph with regions, calibrated image, TCD map, label map, match residual, and image colours against the simulated locus;
- `<out>_summary.csv`: fractions per class or label bin, unassigned fraction, median residuals, and the share of covered pixels whose best alternative fits within 2 ΔE;
- `<out>_maps.npz`: per-pixel `tcd`, `value`, `residual`, `alt_value`, `alt_residual`, and `classes` for whole-number labels;
- `<out>_run.json`: regions, correction matrix and settings.

### `gamma`: measure the camera's decoding exponent

```bash
python python/run_tcd.py gamma --images t5.tif t10.tif t20.tif t40.tif --exposures 5 10 20 40
```

Use one bare-substrate field and change only the exposure time, with auto-exposure, gain and white balance off.

## Using the package

```python
import sys; sys.path.insert(0, "python")
from tcd.system import load_system
from tcd.reference import build_reference, Reference
from tcd.image_tcd import analyse, layer_classes, summary_rows

system = load_system("systems/template.jsonc", {"oxide_nm": "285"})   # constants can be overridden
print(system.n_rows, system.describe())
N, d, names = system.stack(0)                 # indices (471 wavelengths) and thicknesses of candidate 0
R = system.at([{"film_nm": 150, "rough_nm": 0}]).spectra()   # spectrum of any structure

ref = build_reference(system)                 # TMM -> XYZ -> CAM02-UCS for every candidate
ref.save("refs/template_285_D65.csv")
res = analyse("my_image.png", ref, substrate_roi=(10, 10, 40, 40))
print(summary_rows(res))
```

| Module | Content |
| --- | --- |
| `tcd/system.py` | reading and checking system files; candidate rows; layer stacks; materials |
| `tcd/expr.py` | the expression language of system files (parsed, never `eval`-ed) |
| `tcd/tmm.py` | transfer-matrix reflectance: normal, oblique, NA cone average |
| `tcd/materials.py` | n,k library and files, Cauchy, Sellmeier, Maxwell-Garnett, Bruggeman, linear mixing |
| `tcd/colorimetry.py` | spectra to XYZ, CAM02-UCS, CIELAB, sRGB; camera linearisation |
| `tcd/reference.py` | `Reference` (save / load / find / subset), `build_reference`, `ambiguity` |
| `tcd/image_tcd.py` | calibration, per-pixel matching, alternatives, classes, summaries |
| `tcd/plots.py` | reference sheet, spectra figure, result figure |
| `tcd/roi.py`, `tcd/camera.py` | region selection; gamma measurement |
