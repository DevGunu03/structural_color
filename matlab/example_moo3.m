%EXAMPLE_MOO3 Thickness map of exfoliated MoO3 flakes on 100 nm SiO2 / Si.
%   Run from the matlab/ folder. The substrate is taken as the dominant background colour;
%   use 'Substrate', 'click' if the flakes cover most of the field.
here = fileparts(mfilename('fullpath'));
ref = tcd_load_reference(fullfile(here, '..', 'refs', 'moo3_D65.csv'));

[f, p] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Images'}, 'Select a MoO3 micrograph');
if isequal(f, 0), return; end
[~, name] = fileparts(f);
res = tcd_map_image(fullfile(p, f), ref, ...
    'Substrate', 'auto', ...
    'Gamma', 1, ...
    'TMax', 600, ...            % lower it if AFM shows the flakes are thinner: colours repeat every ~130 nm
    'Out', fullfile(here, '..', 'results', name));
disp(res.summary)
% res.altValue / res.altResidual: the best thickness > 40 nm away from the chosen one and how
% well it fits; where altResidual is close to res.residual the colour alone cannot decide.
