%RUN_EXAMPLES Map the bundled example micrographs with the MATLAB code.
%   Run from the matlab/ folder. Results go to ../results/examples/ (figures, summaries, maps).
here = fileparts(mfilename('fullpath'));
img = fullfile(here, '..', 'examples', 'images');
out = fullfile(here, '..', 'results', 'examples');
ps = tcd_load_reference(fullfile(here, '..', 'refs', 'ps_D65.csv'));
moo3 = tcd_load_reference(fullfile(here, '..', 'refs', 'moo3_D65.csv'));

% Substrate colours (sRGB, 0-1) of the two PS microscopes, read off bare-substrate areas;
% the MoO3 images use the dominant background colour instead.
cases = {'ps_chennai_metco', ps,   [0.157 0.263 0.459];
         'ps_olympus',       ps,   [0.580 0.561 0.333];
         'moo3_region_1',    moo3, 'auto';
         'moo3_region_3',    moo3, 'auto';
         'moo3_50x2',        moo3, 'auto'};
for i = 1:size(cases, 1)
    res = tcd_map_image(fullfile(img, [cases{i, 1} '.png']), cases{i, 2}, 'Substrate', cases{i, 3}, ...
        'ShowFigure', false, 'Out', fullfile(out, cases{i, 1}));
    fprintf('\n%s\n', cases{i, 1});
    disp(res.summary(res.summary.fraction > 0.005, :));
end
