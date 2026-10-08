# Python: `run_tcd.py` and the `tcd` package

You need Python 3.10 or newer. Do all commands in the repository folder.

1. Type this command to install the necessary packages:

   ```bash
   pip install -r python/requirements.txt
   ```

2. Type this command to test the installation:

   ```bash
   python python/tests/test_pipeline.py
   ```

   The last line must be `all checks passed`.

To see all options of a command, type `python python/run_tcd.py <command> --help`.

## Commands

| Command | What it does |
| --- | --- |
| `materials` | Shows the materials in the library and the example systems. |
| `add-material` | Adds the optical constants of a material to the library. |
| `spectrum` | Shows the reflectance and the colour of some structures of a system file. |
| `build-ref` | Makes a reference (all structures and their colours) from a system file. |
| `map` | Makes a thickness or layer map of a micrograph, with a reliability score. |
| `gamma` | Measures how the camera encodes its values. |

### `materials`

```bash
python python/run_tcd.py materials                 # the library and the example systems
python python/run_tcd.py materials --show aMoO3    # n and k from 400 nm to 800 nm
```

### `add-material`

```bash
python python/run_tcd.py add-material MoS2 MoS2_downloaded.csv --source "Author et al. 2020, refractiveindex.info"
```

| Option | Use |
| --- | --- |
| `NAME FILE` | The name of the new material and the file with its data. |
| `--source` | Where the data come from. You must give it. |
| `--reference`, `--notes` | A DOI or URL, and other information. |
| `--units` | `auto` (default), `nm` or `um`. |
| `--replace` | Replace an added material with the same name. |

The command saves the material in `data/materials/<NAME>.csv`. For the procedure and the accepted files, refer to [data/materials/README.md](../data/materials/README.md).

### `spectrum`

```bash
python python/run_tcd.py spectrum --system ps --at layers=1,packing=1 --at layers=2,packing=0.5
python python/run_tcd.py spectrum --system systems/my_sample.jsonc --at t=100 --csv
```

The command shows each layer stack, with its colour values. It also makes a figure of the reflectance spectra (`results/<name>_spectra.png`). Give each sweep parameter in `--at`. Without `--at`, the command shows the first, the middle and the last structure.

### `build-ref`

```bash
python python/run_tcd.py build-ref --system moo3                              # makes refs/moo3_D65.csv
python python/run_tcd.py build-ref --system ps --set bead_nm=500 --set oxide_nm=285
python python/run_tcd.py build-ref --system systems/my_sample.jsonc --na 0.5 --spectra
```

| Option | Use |
| --- | --- |
| `--system` | A system file, or the name of a file in `systems/`. |
| `--set NAME=VALUE` | Change a constant of the system file. You can use it more than one time. |
| `--illuminant` | Change the light source: `D65`, `A`, `D50` or a CSV file of your lamp spectrum. |
| `--na` | Change the numerical aperture of the objective. 0 is normal incidence. |
| `--out` | The name of the reference file. The default is `refs/<name>[_<changes>]_<illuminant>.csv`. |
| `--spectra` | Also save all reflectance spectra. |

The command makes three files:

- `<out>.csv`: one row for each structure, with its parameters and its colour values.
- `<out>.json`: the description of the label, the layer stack and the system file.
- `<out>.png`: the reference sheet. Refer to [systems/README.md](../systems/README.md#4-build-the-reference-and-read-its-sheet).

### `map`

```bash
python python/run_tcd.py map --ref refs/ps_D65.csv --interactive
python python/run_tcd.py map --ref refs/moo3_D65.csv --image flakes.png --substrate auto --max-value 350
python python/run_tcd.py map --ref refs/moo3_D65.csv --image flakes.png --substrate auto --check 257nm@249,475,10,10
```

Regions are `x,y,width,height` in pixels. Python counts pixels from 0.

| Option | Use |
| --- | --- |
| `--ref` | The reference file from `build-ref`. |
| `--image` | The micrograph. Without this option, a window opens and you select the file. |
| `--out` | The name of the result files. The default is `results/<image name>`. |
| `--substrate` | The bare-substrate region: `x,y,w,h`, `auto` (the most frequent colour) or `colour:r,g,b`. |
| `--substrate-colour r,g,b` | The colour of the substrate (0–1). The software finds the region. |
| `--interactive` | Select regions with the mouse if you do not give them. |
| `--anchor VALUE@REGION` | A region of known structure. The calibration uses it. Example: `1L@120,40,30,30`. |
| `--check VALUE@REGION` | A region of known thickness (for example from AFM). The software compares it with the map. The calibration does not use it. |
| `--gamma` | How to decode the camera values: a number g, or `srgb`. The default 1 is a linear camera. |
| `--max-value` | Use only structures with a label of this value or less. Example: a maximum thickness from AFM. |
| `--colour-error` | The colour error of the model, in ΔE, for the reliability score. The default comes from the checks, or is 3. |
| `--tolerance` | How near to the correct value a result must be to count as correct. The default is half the gap, or the same layer number. |
| `--min-reliability` | The minimum score for the "only where reliable" panel. The default is 0.5. |
| `--gap` | Two structures that are farther apart than this are different structures. The default comes from the reference. |
| `--max-residual` | Pixels with a larger colour error get no value. The default is 15 ΔE. |
| `--correction` | `auto`, `diagonal` or `matrix`: the type of colour correction. |
| `--sigma`, `--max-side` | The smoothing (default 1.2 px) and the maximum image size (default 1600 px; 0 = no change). |

The command makes these files:

- `<out>.png`: a figure with eight panels. Refer to [Read the results](../README.md#read-the-results).
- `<out>_summary.csv`: the area of each layer number or thickness range, and the share of reliable pixels.
- `<out>_maps.npz`: the values of each pixel: `value`, `reliability`, `fit`, `uniqueness`, `consistency`, `residual`, `alt_value`, `alt_residual`, `tcd`, and `classes` for layer numbers.
- `<out>_run.json`: the regions, the colour correction, the colour error, the checks and all settings.
- `<out>_checks.csv`: the result of each check, if you gave checks.

### `gamma`

```bash
python python/run_tcd.py gamma --images t5.tif t10.tif t20.tif t40.tif --exposures 5 10 20 40
```

1. Put a bare substrate under the microscope.
2. Set the camera to manual exposure, manual gain and manual white balance.
3. Take four to six images. Change only the exposure time.
4. Give the images and the exposure times to the command.

## Use the package in your own scripts

```python
import sys; sys.path.insert(0, "python")
from tcd.system import load_system
from tcd.reference import build_reference
from tcd.image_tcd import analyse, summary_rows

system = load_system("systems/template.jsonc", {"oxide_nm": "285"})
ref = build_reference(system)
res = analyse("my_image.png", ref, substrate_roi=(10, 10, 40, 40),
              checks=[((249, 475, 10, 10), 257.0, "257nm")])
print(summary_rows(res))
print(res.checks)            # how the map compares with the check
print(res.reliability)       # score of each pixel, 0-1
```

| Module | Content |
| --- | --- |
| `tcd/system.py` | Reads and examines system files. Makes the layer stack of each structure. |
| `tcd/expr.py` | The expressions of system files. They are parsed, not executed. |
| `tcd/tmm.py` | The transfer-matrix calculation: normal incidence, oblique incidence, objective NA. |
| `tcd/materials.py` | The n,k library, added materials, files, Cauchy, Sellmeier and mixtures. |
| `tcd/colorimetry.py` | From spectra to XYZ, CAM02-UCS, CIELAB and sRGB. |
| `tcd/reference.py` | References: make, save, load, find a row. |
| `tcd/image_tcd.py` | Calibration, the match of each pixel, layer numbers, summaries. |
| `tcd/reliability.py` | The reliability score and the checks. |
| `tcd/plots.py` | The reference sheet, the spectrum figure and the result figure. |
| `tcd/roi.py`, `tcd/camera.py` | Selection of regions. Measurement of the camera gamma. |
