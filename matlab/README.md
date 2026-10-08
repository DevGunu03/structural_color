# MATLAB: `tcd_*` functions

You need MATLAB R2021a or newer. You do not need toolboxes. Go to this folder before you use the functions. For all options of a function, type `help <function>`.

```matlab
addpath tests; run_tests       % the last line must be "all checks passed"
run_examples                   % maps the five example micrographs
example_own_system             % the full procedure for your own sample, step by step
```

## The procedure

```matlab
tcd_materials                                            % the library and the example systems
sys = tcd_load_system('../systems/my_sample.jsonc');     % read and examine the system file
tcd_spectrum(sys, struct('t', [0 100 200]))              % look at some structures
ref = tcd_build_reference(sys, 'Out', 'auto');           % makes ../refs/my_sample_D65.csv and the sheet
res = tcd_map_image('', ref, 'Out', '../results/my_map'); % select the image, click the substrate
disp(res.summary)
disp(res.checks)                                         % if you gave 'Checks'
```

For the system files, refer to [../systems/README.md](../systems/README.md). The example systems load by name, for example `tcd_build_reference('moo3')`.

## Functions

| Function | What it does |
| --- | --- |
| `tcd_load_system` | Reads a system file. Examines all keys, materials and expressions. Makes the list of structures. `'Set'` changes constants. |
| `tcd_spectrum` | Shows the reflectance, the colour and the layer stack of some structures. |
| `tcd_build_reference` | Calculates the colour of each structure. `'Out'` saves the reference and the reference sheet. |
| `tcd_plot_reference` | Draws the reference sheet. |
| `tcd_map_image` | Makes a thickness or layer map of a micrograph, with a reliability score and optional checks. |
| `tcd_add_material` | Adds the optical constants of a material to the library (`data/materials/`). |
| `tcd_materials` | Shows the library. Gives n + ik of a material. |
| `tcd_tmm` | The transfer-matrix calculation. |
| `tcd_load_reference`, `tcd_save_reference` | Read and write references. Python uses the same files. |
| `tcd_measure_gamma` | Measures the camera gamma from images with different exposure times. |
| `tcd_colour` | Colour conversions. |

The folder `private/` contains the helper functions.

## Options of `tcd_map_image`

Regions are `[x y width height]` in pixels. MATLAB counts pixels from 1.

| Option | Use |
| --- | --- |
| `'Substrate'` | The bare-substrate region: `[x y w h]`, `'click'` (default), `'auto'` or a colour `[r g b]`. |
| `'Anchors'` | Regions of known structure for the calibration: `{'1L', [x y w h]; ...}`. |
| `'Checks'` | Regions of known thickness for comparison only: `{'257nm', [x y w h]; ...}`. |
| `'Gamma'` | How to decode the camera values: a number, or `'srgb'`. The default 1 is a linear camera. |
| `'MaxValue'` | Use only structures with a label of this value or less. |
| `'ColourError'` | The colour error of the model in ΔE. The default comes from the checks, or is 3. |
| `'Tolerance'` | How near to the correct value a result must be to count as correct. |
| `'MinReliability'` | The minimum score for the "only where reliable" panel. The default is 0.5. |
| `'Out'` | The name of the result files. Without it, the function saves no files. |

## Results of `tcd_map_image`

| Field | Content |
| --- | --- |
| `value` | The label of each pixel (thickness, layer number, ...). NaN means no match. |
| `reliability` | The score of each pixel, from 0 to 1: `fit .* uniqueness .* consistency`. |
| `fit`, `uniqueness`, `consistency` | The three parts of the score. |
| `checks` | A table: the result of each check. |
| `sigma` | The colour error that the score uses: camera noise, model error, total. |
| `classes` | Layer numbers only: the layer number of each pixel. −1 means no match. |
| `residual` | The colour difference between each pixel and its match. |
| `altValue`, `altResidual` | The best different structure, and its colour difference. |
| `tcd` | The colour difference from the substrate. |
| `summary` | A table of the areas and the share of reliable pixels. |

> **CAUTION:** A high reliability score does not prove that a thickness is correct. If the model is wrong (for example an incorrect oxide thickness), a wrong thickness can get a high score. Use `'Checks'` to find these errors.

## The transfer-matrix code

`tcd_tmm.m` replaces `legacy/TransferMatrix/TransferMatrix_packing.m` and `TransferMatrix_multiple.m`. It uses the same matrices. The tests compare it with spectra from the old code. The changes are:

- It calculates all wavelengths at the same time.
- A system file gives the layers and the mixtures. The old code used the Maxwell-Garnett equation only for layers with the name `'PS-beads'`.
- It can calculate oblique incidence and the light cone of an objective (`'NA'`).
- It gives the true reflectance. The old code gave R·D65 / (12·max(R·D65)), which changed the scale of each structure differently.

The header of `tcd_tmm.m` explains the method and the changes in more detail.
