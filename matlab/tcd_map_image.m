function res = tcd_map_image(imageFile, ref, opts)
%TCD_MAP_IMAGE Layer-number / thickness map of an optical micrograph from a simulated reference.
%
%   res = tcd_map_image('flakes.png', tcd_load_reference('../refs/moo3_D65.csv'))   % click the substrate
%   res = tcd_map_image('beads.png', '../refs/ps_D65.csv', 'Substrate', [187 434 40 40], 'Out', 'results/beads')
%   res = tcd_map_image('flakes.png', ref, 'Substrate', 'auto', 'TMax', 350)
%
%   Steps: read + (box) downsample + Gaussian smoothing -> linearise camera values (Gamma) ->
%   calibrate on the bare substrate (and optional anchors of known structure) so it reproduces the
%   simulated substrate colour -> CAM02-UCS per pixel -> TCD map (Delta E from substrate, absolute)
%   -> nearest reference row in J'a'b' for every pixel.
%
%   Regions are [x y w h] in pixels of the original image, x and y = 1-based top-left corner.
%
%   Options (name-value):
%     Substrate ('click')   [x y w h] | 'click' | 'auto' (dominant background colour) | [r g b] in 0-1
%     Anchors ({})          N-by-2 cell of {label, region}, label '1L', '2L' (PS) or '250nm' (MoO3);
%                           only for regions whose structure is known (SEM/AFM)
%     Gamma (1)             camera decoding: number g (intensity = value^g) or 'srgb'
%     Sigma (1.2)           Gaussian smoothing in pixels
%     MaxResidual (15)      Delta E beyond which a pixel is left unassigned
%     MaxSide (1600)        downsample longer side to this many pixels (0 = never)
%     Correction ('auto')   'auto' | 'diagonal' | 'matrix'
%     Gap (40)              MoO3: report the best alternative more than Gap nm away
%     TMax ([])             MoO3: only consider thicknesses <= TMax (prior from AFM)
%     MaxLayers ([])        PS: only consider up to this many layers
%     Out ('')              path stem: writes <Out>.png, <Out>_summary.csv, <Out>_maps.mat
%     ShowFigure (true)
%
%   res fields: tcd, value (nm or effective layer number), residual, index (0 = unassigned),
%   layers (PS, majority-filtered, -1 = unassigned), altValue/altResidual (MoO3), calibrated,
%   correction, regions, summary.
%
%   See also TCD_BUILD_REFERENCE, TCD_LOAD_REFERENCE, TCD_MEASURE_GAMMA.
arguments
    imageFile (1,:) char
    ref
    opts.Substrate = 'click'
    opts.Anchors cell = {}
    opts.Gamma = 1
    opts.Sigma (1,1) double = 1.2
    opts.MaxResidual (1,1) double = 15
    opts.MaxSide (1,1) double = 1600
    opts.Correction (1,:) char {mustBeMember(opts.Correction, {'auto', 'diagonal', 'matrix'})} = 'auto'
    opts.Gap (1,1) double = 40
    opts.TMax = []
    opts.MaxLayers = []
    opts.Out (1,:) char = ''
    opts.ShowFigure (1,1) logical = true
end
if ischar(ref) || isstring(ref), ref = tcd_load_reference(char(ref)); end
if strcmp(ref.system, 'MoO3') && ~isempty(opts.TMax), ref = subset(ref, ref.value <= opts.TMax); end
if strcmp(ref.system, 'PS') && ~isempty(opts.MaxLayers), ref = subset(ref, ref.params.layers <= opts.MaxLayers); end
c = colour_const();

[img, scale] = load_image(imageFile, opts.MaxSide);
[H, W, ~] = size(img);
saturated = any(img >= 0.985, 3);
work = gauss_smooth(img, opts.Sigma);
if ischar(opts.Gamma) || isstring(opts.Gamma)
    lin = srgb_codec(work, 'decode');
else
    lin = work .^ opts.Gamma;
end

% calibration regions
regions = struct('name', 'substrate', 'roi', resolve_roi(opts.Substrate, img, scale, 'bare-substrate', 40), 'row', 1);
for a = 1:size(opts.Anchors, 1)
    regions(end + 1) = struct('name', opts.Anchors{a, 1}, ...
        'roi', resolve_roi(opts.Anchors{a, 2}, img, scale, opts.Anchors{a, 1}, 12), ...
        'row', find_row(ref, opts.Anchors{a, 1})); %#ok<AGROW>
end
measured = zeros(numel(regions), 3);
for k = 1:numel(regions)
    measured(k, :) = median(roi_pixels(lin, regions(k).roi, scale), 1);
end
target = (ref.XYZ([regions.row], :) / 100) * c.XYZ2RGB';
M = fit_correction(measured, target, opts.Correction);

% colour of every pixel, TCD map and assignment
linCal = max(reshape(lin, [], 3) * M', 0);
Jab = xyz_to_cam02ucs(100 * linCal * c.RGB2XYZ');
Jab(~isfinite(Jab)) = 0;
tcd = vecnorm(Jab - ref.Jab(1, :), 2, 2);
if strcmp(ref.system, 'MoO3')
    [idx, dist, altVal, altDist] = nearest_ref(Jab, ref, ~saturated(:), opts.Gap);
else
    [idx, dist] = nearest_ref(Jab, ref);
end
bad = dist > opts.MaxResidual | saturated(:);
idx(bad) = 0;
value = nan(H * W, 1);
value(~bad) = ref.value(idx(~bad));

res.tcd = reshape(tcd, H, W);
res.value = reshape(value, H, W);
res.residual = reshape(dist, H, W);
res.index = reshape(idx, H, W);
res.saturated = saturated;
if strcmp(ref.system, 'MoO3')
    altVal(bad) = nan;  altDist(bad) = nan;
    res.altValue = reshape(altVal, H, W);
    res.altResidual = reshape(altDist, H, W);
else
    nCls = floor(max(ref.value)) + 1;
    cls = floor(res.value + 0.5);                     % half up, as in the Python code
    cls(isnan(cls)) = -1;
    res.layers = mode_filter(cls, 5, nCls);
end
res.image = img;
res.calibrated = min(max(srgb_codec(min(reshape(linCal, H, W, 3), 1), 'encode'), 0), 1);
res.Jab = reshape(Jab, H, W, 3);
res.correction = M;
res.regions = regions;
res.scale = scale;
res.reference = ref;
res.settings = rmfield(opts, 'ShowFigure');
res.settings.image = imageFile;
res.summary = summarise(res, ref);

if opts.ShowFigure || ~isempty(opts.Out)
    fig = plot_result(res, opts.ShowFigure);
    if ~isempty(opts.Out)
        folder = fileparts(opts.Out);
        if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
        exportgraphics(fig, [opts.Out '.png'], 'Resolution', 150);
        writetable(res.summary, [opts.Out '_summary.csv']);
        maps = rmfield(res, {'image', 'calibrated', 'Jab', 'reference', 'saturated'}); %#ok<NASGU>
        save([opts.Out '_maps.mat'], '-struct', 'maps');
    end
    if ~opts.ShowFigure, close(fig); end
end
end

% -------------------------------------------------------------------------------------------
function px = roi_pixels(img, roi, scale)
x0 = round((roi(1) - 1) * scale) + 1;  y0 = round((roi(2) - 1) * scale) + 1;
w = max(round(roi(3) * scale), 1);       h = max(round(roi(4) * scale), 1);
rows = max(y0, 1):min(y0 + h - 1, size(img, 1));
cols = max(x0, 1):min(x0 + w - 1, size(img, 2));
if isempty(rows) || isempty(cols), error('tcd:roi', 'Region [%s] lies outside the image', num2str(roi)); end
px = reshape(img(rows, cols, :), [], 3);
end

function row = find_row(ref, label)
s = lower(strrep(label, ' ', ''));
if any(strcmp(s, {'substrate', '0', '0l', '0nm'}))
    row = 1;
elseif strcmp(ref.system, 'PS') && endsWith(s, 'l')
    n = str2double(s(1:end - 1));
    row = find(ref.params.layers == n & abs(ref.params.packing - 1) < 1e-9, 1);
else
    [~, row] = min(abs(ref.value - str2double(erase(s, 'nm'))));
end
if isempty(row), error('tcd:anchor', 'No reference row for anchor "%s"', label); end
end

function ref = subset(ref, keep)
keep(1) = true;                                   % always keep the substrate
f = fieldnames(ref.params);
for i = 1:numel(f), ref.params.(f{i}) = ref.params.(f{i})(keep); end
for g = {'value', 'XYZ', 'Jab', 'Lab', 'sRGB', 'dE'}
    ref.(g{1}) = ref.(g{1})(keep, :);
end
if isfield(ref, 'spectra'), ref.spectra = ref.spectra(keep, :); end
end

function T = summarise(res, ref)
total = numel(res.index);
items = {'unassigned'};  px = sum(res.index(:) == 0);  med = nan;
if strcmp(ref.system, 'PS')
    cls = floor(res.value + 0.5);
    for k = 0:floor(max(ref.value))
        m = cls == k;
        if k == 0, items{end + 1} = 'substrate'; else, items{end + 1} = sprintf('%dL', k); end %#ok<AGROW>
        px(end + 1) = sum(m(:)); med(end + 1) = median(res.residual(m)); %#ok<AGROW>
    end
else
    flake = res.index > 0 & res.value > 20;
    close2 = flake & (res.altResidual - res.residual < 2);
    items{end + 1} = 'flake px with alternative thickness within 2 dE';
    px(end + 1) = sum(close2(:));  med(end + 1) = nan;
    edges = 0:50:max(ref.value) + 50;
    for k = 1:numel(edges) - 1
        m = res.value >= edges(k) & res.value < edges(k + 1);
        items{end + 1} = sprintf('%d-%d nm', edges(k), edges(k + 1)); %#ok<AGROW>
        px(end + 1) = sum(m(:)); med(end + 1) = median(res.residual(m)); %#ok<AGROW>
    end
end
frac = px(:) / total;
if strcmp(ref.system, 'MoO3'), frac(2) = px(2) / max(sum(flake(:)), 1); end
T = table(items(:), px(:), frac, med(:), 'VariableNames', {'item', 'pixels', 'fraction', 'median_residual'});
end
