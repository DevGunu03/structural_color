%EXAMPLE_PS Layer-number map of a PS-bead micrograph (300 nm beads on 100 nm SiO2 / Si).
%   Run from the matlab/ folder. Pick an image, then click two opposite corners of a
%   bare-substrate region. Results go to ../results/<image name>.*
here = fileparts(mfilename('fullpath'));
ref = tcd_load_reference(fullfile(here, '..', 'refs', 'ps_D65.csv'));
% Other beads or oxide: build a reference from the same system file with other constants, e.g.
%   ref = tcd_build_reference('ps', 'Set', {'bead_nm', 500; 'oxide_nm', 285}, 'Out', 'auto');
% or try the close-packed model:  ref = tcd_build_reference('ps_hcp', 'Out', 'auto');

[f, p] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Images'}, 'Select a PS-bead micrograph');
if isequal(f, 0), return; end
[~, name] = fileparts(f);
res = tcd_map_image(fullfile(p, f), ref, ...
    'Substrate', 'click', ...   % or [x y w h], or the substrate colour [r g b]
    'Gamma', 1, ...             % measure it once per camera with tcd_measure_gamma
    'Out', fullfile(here, '..', 'results', name));
disp(res.summary)
% res.classes: layer number per pixel (-1 = unassigned); res.value: effective layer number.
