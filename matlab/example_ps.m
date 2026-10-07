%EXAMPLE_PS Layer-number map of a PS-bead micrograph (300 nm beads on 100 nm SiO2 / Si).
%   Run from the matlab/ folder. Pick an image, then click two opposite corners of a
%   bare-substrate region. Results go to ../results/<image name>.*
here = fileparts(mfilename('fullpath'));
ref = tcd_load_reference(fullfile(here, '..', 'refs', 'ps_D65.csv'));
% For another stack, build your own reference instead, e.g.
%   ref = tcd_build_reference('PS', 'Bead', 500, 'Oxide', 285, 'Out', fullfile(here, '..', 'refs', 'ps_500nm_285ox.csv'));

[f, p] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Images'}, 'Select a PS-bead micrograph');
if isequal(f, 0), return; end
[~, name] = fileparts(f);
res = tcd_map_image(fullfile(p, f), ref, ...
    'Substrate', 'click', ...   % or [x y w h], or the substrate colour [r g b]
    'Gamma', 1, ...             % measure it once per camera with tcd_measure_gamma
    'Out', fullfile(here, '..', 'results', name));
disp(res.summary)
