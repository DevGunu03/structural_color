%EXAMPLE_OWN_SYSTEM The whole workflow for a sample of your own, step by step.
%   Run from the matlab/ folder, one section at a time (Ctrl+Enter). The example uses the
%   bundled template (a polymer film on SiO2 / Si); replace SYSTEM with your own file.
%   Guide to writing system files: ../systems/README.md

%% 1. Describe the sample in a system file
% Copy ../systems/template.jsonc to ../systems/my_sample.jsonc and edit it: the layers from the
% light side down, their materials, the substrate, and the parameters that vary (the sweep).
% Materials available by name:
tcd_materials
here = fileparts(mfilename('fullpath'));
SYSTEM = fullfile(here, '..', 'systems', 'template.jsonc');

%% 2. Load it: this checks every key, material and expression
sys = tcd_load_system(SYSTEM);                  % constants can be changed: 'Set', {'oxide_nm', 285}
fprintf('%s: %d candidate structures, label "%s" (%s)\n', sys.name, sys.nRows, sys.label.name, sys.label.unit);
disp(sys.rows)

%% 3. Look at a few single structures before building everything
% Are the layers in the right order, the thicknesses what you meant, the colours plausible?
tcd_spectrum(sys, struct('film_nm', [0 100 150 150], 'rough_nm', [0 0 0 20]))

%% 4. Build the colour reference (one TMM spectrum + colour per candidate) and its sheet
out = fullfile(here, '..', 'refs', [sys.name '_D65.csv']);
ref = tcd_build_reference(sys, 'Out', out);     % writes .csv, .json and the reference sheet .png
tcd_plot_reference(ref);                        % where colours repeat, the map cannot decide
% If your objective has a large NA, compare with a cone-averaged reference:
%   ref = tcd_build_reference(sys, 'NA', 0.5, 'Out', 'auto');

%% 5. Map a micrograph with it
% Choose the image, then click two corners of a bare-substrate region (the structure of row 1).
% Measure your camera's gamma once with tcd_measure_gamma; 1 means it writes linear data.
res = tcd_map_image('', ref, 'Substrate', 'click', 'Gamma', 1, ...
    'Out', fullfile(here, '..', 'results', [sys.name '_map']));
disp(res.summary)
% res.value     label of every pixel (NaN = unassigned)
% res.residual  how far the pixel colour is from the matched simulated colour (Delta E)
% res.altValue  the best match more than label.gap away; if res.altResidual is close to
%               res.residual, the colour alone cannot tell the two apart
