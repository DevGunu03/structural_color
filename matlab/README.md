# MATLAB: `tcd_*` functions

MATLAB R2021a or newer (tested on R2025b). No toolboxes, and no CIECAM02 download. Work from this folder. `help <function>` documents every option.

```matlab
addpath tests; run_tests       % prints "all checks passed"
run_examples                   % maps the five bundled micrographs
example_own_system             % the whole workflow for your own sample, section by section
```

## The workflow

```matlab
tcd_materials                                            % n,k library and bundled systems
sys = tcd_load_system('../systems/my_sample.jsonc');     % read + check the system file
tcd_spectrum(sys, struct('t', [0 100 200]))              % check single structures
ref = tcd_build_reference(sys, 'Out', 'auto');           % ../refs/my_sample_D65.csv + .json + sheet .png
res = tcd_map_image('', ref);                            % choose the image, click a bare-substrate region
disp(res.summary)
```

System files are described in [../systems/README.md](../systems/README.md). The bundled ones load by name: `tcd_build_reference('moo3')`, `tcd_build_reference('ps', 'Set', {'oxide_nm', 285})`.

## Functions

| Function | Does |
| --- | --- |
| `tcd_load_system` | reads a system file (`//` comments allowed) and checks every key, material and expression; expands the sweep into candidate rows. `'Set'` changes constants |
| `tcd_spectrum` | reflectance spectrum, colour and stack of chosen structures; plots them when called without output |
| `tcd_build_reference` | TMM spectrum and colour (XYZ, CAM02-UCS, CIELAB, sRGB) of every candidate; `'Out'` saves CSV + JSON + reference sheet. `'Illuminant'`, `'NA'` override the file |
| `tcd_plot_reference` | the reference sheet: colour strip, ΔE from the substrate, ΔE to the nearest look-alike, a′b′ path, colour chart |
| `tcd_map_image` | calibrates a micrograph on the bare substrate (and optional anchors) and assigns every pixel to the nearest candidate; works with the reference of any system |
| `tcd_tmm` | the transfer-matrix reflectance itself: normal incidence, one oblique angle, or the cone of an objective (`'NA'`) |
| `tcd_materials` | lists the library; returns n + ik of a library name, an n,k file or a constant |
| `tcd_load_reference`, `tcd_save_reference` | references in the CSV + JSON format shared with Python |
| `tcd_measure_gamma` | camera decoding exponent from an exposure series |
| `tcd_colour` | colour conversions (sRGB, XYZ, CAM02-UCS, CIELAB, reflectance to XYZ) |

`private/` holds the helpers:

- `expr_eval`: the expression language;
- `system_stack`: one candidate's layer stack;
- `read_jsonc`, `load_nk` and `ema_mix`: reading files, n,k data and mixtures;
- the colour conversions;
- image I/O and region selection;
- the nearest-reference search and the figures.

## The transfer-matrix code

`tcd_tmm.m` replaces `legacy/TransferMatrix/TransferMatrix_packing.m` and `TransferMatrix_multiple.m`. The matrices and the reflectance are the same, checked against spectra saved from the old code. The differences:

- **Vectorised over wavelength.**
- **Any stack.** Materials, thicknesses and mixtures come from a system file. The old code applied Maxwell-Garnett only to layers named `'PS-beads'` and read the fill fractions with `eval('vol_incl_' + i)`.
- **Oblique incidence and the NA cone average.**
- **R is the true reflectance.** The lamp and the observer are applied once, in the colorimetry. The old code returned R·D65 / (12·max(R·D65)), so every structure was rescaled differently.

The header of `tcd_tmm.m` explains the method and these changes in full. The field-profile, absorption and photocurrent parts of the original are not needed for colour and remain in `legacy/TransferMatrix/`.

## `tcd_map_image` results

| Field | Content |
| --- | --- |
| `value` | label of the matched candidate for every pixel (thickness, effective layer number, ...); NaN = unassigned |
| `classes` | whole-number labels only: classes after rounding half up and a 5×5 majority filter; −1 = unassigned |
| `residual` | ΔE between the pixel and its match |
| `altValue`, `altResidual` | best candidate more than `label.gap` away, and its ΔE; close to `residual` means the colour cannot decide |
| `tcd` | ΔE from the simulated substrate (absolute) |
| `summary` | table of fractions per class or label bin, unassigned share, ambiguous share |
| `calibrated`, `correction`, `regions`, `settings` | calibrated image, colour-correction matrix, regions used, options |

Regions are `[x y w h]` with a **1-based** top-left corner (Python uses 0-based).
