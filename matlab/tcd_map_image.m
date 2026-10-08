function res = tcd_map_image(imageFile, ref, opts)
%TCD_MAP_IMAGE Thickness / layer-number map of an optical micrograph from a simulated reference.
%
%   res = tcd_map_image('', '../refs/moo3_D65.csv')                 % pick the image, click the substrate
%   res = tcd_map_image('beads.png', '../refs/ps_D65.csv', 'Substrate', [187 434 40 40], 'Out', 'results/beads')
%   res = tcd_map_image('flakes.png', ref, 'Substrate', 'auto', 'MaxValue', 350)
%
%   Works with the reference of any system file (tcd_build_reference): what the map reports
%   (thickness, layer number, ...) and how it is shown come from the reference's label.
%
%   Steps: read + (box) downsample + Gaussian smoothing -> linearise camera values (Gamma) ->
%   calibrate on the bare substrate (and optional anchors of known structure) so it reproduces the
%   simulated colour of reference row 1 -> CAM02-UCS per pixel -> TCD map (Delta E from the
%   substrate, absolute) -> nearest reference row in J'a'b' for every pixel, plus the best
%   alternative more than Gap away (colours repeat with thickness) -> for whole-number labels
%   (label.classes), classes rounded half up and cleaned with a 5x5 majority filter ->
%   reliability score per pixel (fit x uniqueness x consistency, 0-1) and optional checks
%   against regions of known thickness.
%
%   Reliability: fit = how well the matched structure explains the colour; uniqueness = how much
%   of the evidence points within Tolerance of the reported value rather than to a look-alike;
%   consistency = how many neighbours agree. The colour error it assumes is the camera noise
%   (measured in the substrate region) plus the model error (ColourError; by default estimated
%   from the Checks, or 3 Delta E). The score cannot see an error of the model itself (wrong
%   oxide thickness or optical constants), so add Checks measured by AFM whenever you can.
%
%   Regions are [x y w h] in pixels of the original image, x and y = 1-based top-left corner.
%
%   Options (name-value):
%     Substrate ('click')   [x y w h] | 'click' | 'auto' (dominant background colour) | [r g b] in 0-1
%     Anchors ({})          N-by-2 cell of {label, region}; label '1L', '250nm', '2' or
%                           'layers=2,packing=1'; only for regions whose structure is known (SEM/AFM/Raman)
%     Gamma (1)             camera decoding: number g (intensity = value^g) or 'srgb'
%     Sigma (1.2)           Gaussian smoothing in pixels
%     MaxResidual (15)      Delta E beyond which a pixel is left unassigned
%     MaxSide (1600)        downsample longer side to this many pixels (0 = never)
%     Correction ('auto')   'auto' | 'diagonal' | 'matrix'
%     Gap ([])              ambiguity gap in label units (default: the reference's label.gap)
%     Checks ({})           N-by-2 cell of {known value, region}, e.g. {'230nm', [136 439 10 10]}:
%                           regions measured another way (AFM). Compared with the map, never used
%                           for calibration
%     ColourError ([])      model colour error (Delta E per axis) for the reliability score;
%                           default: estimated from the Checks, else 3
%     Tolerance ([])        how close (label units) counts as right; default half the gap, or the class
%     MinReliability (0.5)  pixels below it are greyed out in the 'reliable only' panel
%     MaxValue ([])         only consider candidates whose label is <= MaxValue (prior from AFM, ...)
%                           (TMax and MaxLayers are accepted as older names for it)
%     Out ('')              path stem: writes <Out>.png, <Out>_summary.csv, <Out>_maps.mat
%     ShowFigure (true)
%
%   res fields: tcd, value (label of the nearest reference), residual, index (0 = unassigned),
%   altValue / altResidual (best alternative more than Gap away), classes (label.classes
%   references; -1 = unassigned), reliability, fit, uniqueness, consistency, sigma (colour
%   error used), checks (table), calibrated, correction, regions, summary.
%
%   See also TCD_BUILD_REFERENCE, TCD_LOAD_REFERENCE, TCD_MEASURE_GAMMA.
arguments
    imageFile
    ref
    opts.Substrate = 'click'
    opts.Anchors cell = {}
    opts.Gamma = 1
    opts.Sigma (1,1) double = 1.2
    opts.MaxResidual (1,1) double = 15
    opts.MaxSide (1,1) double = 1600
    opts.Correction (1,:) char {mustBeMember(opts.Correction, {'auto', 'diagonal', 'matrix'})} = 'auto'
    opts.Gap = []
    opts.Checks cell = {}
    opts.ColourError = []
    opts.Tolerance = []
    opts.MinReliability (1,1) double = 0.5
    opts.MaxValue = []
    opts.TMax = []
    opts.MaxLayers = []
    opts.Out (1,:) char = ''
    opts.ShowFigure (1,1) logical = true
end
if ischar(ref) || isstring(ref), ref = tcd_load_reference(char(ref)); end
if isempty(imageFile)
    [f, p] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Images'}, 'Select a micrograph');
    if isequal(f, 0), error('tcd:image', 'No image selected'); end
    imageFile = fullfile(p, f);
end
imageFile = char(imageFile);
maxValue = [opts.MaxValue, opts.TMax, opts.MaxLayers];
if ~isempty(maxValue), ref = subset(ref, ref.value <= maxValue(1)); end
gap = opts.Gap;
if isempty(gap), gap = ref.label.gap; end
tol = opts.Tolerance;
if isempty(tol)
    if ref.label.classes, tol = 0.5; else, tol = ref.label.gap / 2; end
end
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
        'row', ref_find(ref, opts.Anchors{a, 1})); %#ok<AGROW>
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
[idx, dist] = nearest_ref(Jab, ref);
bad = dist > opts.MaxResidual | saturated(:);
value = nan(H * W, 1);
value(~bad) = ref.value(idx(~bad));
valueMap = reshape(value, H, W);
distMap = reshape(dist, H, W);

% colour error for the reliability score: camera noise (substrate) + model error (checks or default)
J3 = reshape(Jab, H, W, 3);
sub = roi_pixels(J3, regions(1).roi, scale);
noise = sqrt(mean(sum((sub - median(sub, 1)) .^ 2, 2)) / 3);
checks = struct('name', {}, 'known', {}, 'roi', {}, 'rows', {}, 'cols', {});
for k = 1:size(opts.Checks, 1)
    roi = resolve_roi(opts.Checks{k, 2}, img, scale, ['check ' opts.Checks{k, 1}], 12);
    [rr, cc] = roi_index(roi, scale, H, W);
    checks(k) = struct('name', opts.Checks{k, 1}, 'known', ref_value(ref, opts.Checks{k, 1}), 'roi', roi, 'rows', rr, 'cols', cc);
end
rightRes = [];
for k = 1:numel(checks)
    v = valueMap(checks(k).rows, checks(k).cols);  r = distMap(checks(k).rows, checks(k).cols);
    ok = ~isnan(v) & same_structure(ref, v, checks(k).known, tol);
    rightRes = [rightRes; r(ok)]; %#ok<AGROW>
end
if ~isempty(opts.ColourError)
    model = opts.ColourError;  source = 'given';
elseif numel(rightRes) >= 30
    model = sqrt(max((median(rightRes) / 1.5381722) ^ 2 - noise ^ 2, 0.25));  source = 'checks';
else
    model = 3;  source = 'default';
end
sigma = hypot(noise, model);

[~, ~, altVal, altDist, uniq] = nearest_ref(Jab, ref, ~bad, gap, tol, sigma);
idx(bad) = 0;
altVal(bad) = nan;  altDist(bad) = nan;  uniq(bad) = nan;
[fit, cons] = reliability_parts(valueMap, distMap, ref, tol, sigma);
uniq = reshape(uniq, H, W);
score = fit .* uniq .* cons;

res.tcd = reshape(tcd, H, W);
res.value = reshape(value, H, W);
res.residual = reshape(dist, H, W);
res.index = reshape(idx, H, W);
res.altValue = reshape(altVal, H, W);
res.altResidual = reshape(altDist, H, W);
res.reliability = score;
res.fit = fit;
res.uniqueness = uniq;
res.consistency = cons;
res.sigma = struct('noise', noise, 'model', model, 'total', sigma, 'model_source', source);
res.checks = check_table(checks, valueMap, distMap, score, ref, tol);
res.saturated = saturated;
if ref.label.classes
    res.classes = classes_of(res.value, ref, 5);
end
res.image = img;
res.calibrated = min(max(srgb_codec(min(reshape(linCal, H, W, 3), 1), 'encode'), 0), 1);
res.Jab = reshape(Jab, H, W, 3);
res.correction = M;
res.regions = regions;
res.scale = scale;
res.reference = ref;
res.settings = rmfield(opts, {'ShowFigure', 'TMax', 'MaxLayers', 'Checks'});
res.settings.Gap = gap;
res.settings.Tolerance = tol;
res.settings.image = imageFile;
res.summary = summarise(res, ref, gap);

if opts.ShowFigure || ~isempty(opts.Out)
    fig = plot_result(res, opts.ShowFigure);
    if ~isempty(opts.Out)
        folder = fileparts(opts.Out);
        if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
        exportgraphics(fig, [opts.Out '.png'], 'Resolution', 150);
        writetable(res.summary, [opts.Out '_summary.csv']);
        if ~isempty(res.checks), writetable(res.checks, [opts.Out '_checks.csv']); end
        maps = rmfield(res, {'image', 'calibrated', 'Jab', 'reference', 'saturated'}); %#ok<NASGU>
        save([opts.Out '_maps.mat'], '-struct', 'maps');
    end
    if ~opts.ShowFigure, close(fig); end
end
for k = 1:height(res.checks)
    c = res.checks(k, :);
    verdict = 'OK';  if c.within_tolerance < 0.5, verdict = 'DISAGREES'; end
    fprintf('check %s: map median %.4g %s vs known %g; %.0f %% of pixels within tolerance, median reliability %.2f  [%s]\n', ...
        c.check{1}, c.median_value, ref.label.unit, c.known, 100 * c.within_tolerance, c.median_reliability, verdict);
end
if ~isempty(res.checks) && all(res.checks.within_tolerance < 0.5)
    warning('tcd:checks', ['The map disagrees with every check. Do not trust it: check the system file ' ...
        '(oxide thickness, optical constants), the camera gamma and the substrate region.']);
end
end

% -------------------------------------------------------------------------------------------
function [rows, cols] = roi_index(roi, scale, H, W)
x0 = round((roi(1) - 1) * scale) + 1;  y0 = round((roi(2) - 1) * scale) + 1;
w = max(round(roi(3) * scale), 1);       h = max(round(roi(4) * scale), 1);
rows = max(y0, 1):min(y0 + h - 1, H);
cols = max(x0, 1):min(x0 + w - 1, W);
if isempty(rows) || isempty(cols), error('tcd:roi', 'Region [%s] lies outside the image', num2str(roi)); end
end

function T = check_table(checks, valueMap, distMap, score, ref, tol)
% How the map compares with each region of known structure.
n = numel(checks);
[med, err, within, rel, res, assigned, px] = deal(nan(n, 1));
for k = 1:n
    v = valueMap(checks(k).rows, checks(k).cols);  r = distMap(checks(k).rows, checks(k).cols);
    s = score(checks(k).rows, checks(k).cols);
    a = ~isnan(v);
    px(k) = numel(v);  assigned(k) = sum(a(:));
    if any(a(:))
        med(k) = median(v(a));  err(k) = median(abs(v(a) - checks(k).known));
        within(k) = mean(same_structure(ref, v(a), checks(k).known, tol));
        rel(k) = median(s(a));  res(k) = median(r(a));
    end
end
T = table({checks.name}', [checks.known]', reshape(cat(1, checks.roi), [], 4), px, assigned, med, err, within, rel, res, ...
    'VariableNames', {'check', 'known', 'roi', 'pixels', 'assigned', 'median_value', 'median_abs_error', ...
    'within_tolerance', 'median_reliability', 'median_residual'});
end

function px = roi_pixels(img, roi, scale)
x0 = round((roi(1) - 1) * scale) + 1;  y0 = round((roi(2) - 1) * scale) + 1;
w = max(round(roi(3) * scale), 1);       h = max(round(roi(4) * scale), 1);
rows = max(y0, 1):min(y0 + h - 1, size(img, 1));
cols = max(x0, 1):min(x0 + w - 1, size(img, 2));
if isempty(rows) || isempty(cols), error('tcd:roi', 'Region [%s] lies outside the image', num2str(roi)); end
px = reshape(img(rows, cols, :), [], 3);
end

function cls = classes_of(value, ref, clean)
% Whole classes (rounded half up, as in the Python code), -1 = unassigned, majority-filtered.
lo = floor(min(ref.value) + 0.5);
n = floor(max(ref.value) + 0.5) - lo + 1;
cls = floor(value + 0.5) - lo;
cls(isnan(cls)) = -1;
if clean > 1, cls = mode_filter(cls, clean, n); end
cls(cls >= 0) = cls(cls >= 0) + lo;
end

function ref = subset(ref, keep)
keep(1) = true;                                   % always keep the calibration row
f = fieldnames(ref.params);
for i = 1:numel(f), ref.params.(f{i}) = ref.params.(f{i})(keep); end
for g = {'value', 'XYZ', 'Jab', 'Lab', 'sRGB', 'dE'}
    ref.(g{1}) = ref.(g{1})(keep, :);
end
if isfield(ref, 'spectra'), ref.spectra = ref.spectra(keep, :); end
end

function T = summarise(res, ref, gap)
% Pixel fractions per class or per label bin, plus how often the colour is ambiguous
% ('covered' = assigned and more than gap/2 from the calibration row's label).
total = numel(res.index);
v = res.value;
unit = ref.label.unit;
items = {'unassigned'};  px = sum(res.index(:) == 0);  frac = px / total;  med = nan;
vv = v; vv(isnan(vv)) = ref.value(1);
covered = res.index > 0 & abs(vv - ref.value(1)) > gap / 2;
close2 = covered & (res.altResidual - res.residual < 2);
items{end + 1} = strtrim(regexprep(sprintf('covered px with an alternative > %g %s away within 2 dE', gap, unit), '\s+', ' '));
px(end + 1) = sum(close2(:));  frac(end + 1) = px(end) / max(sum(covered(:)), 1);  med(end + 1) = nan;
reliable = covered & res.reliability >= 0.5;
items{end + 1} = 'covered px with reliability >= 0.5';
px(end + 1) = sum(reliable(:));  frac(end + 1) = px(end) / max(sum(covered(:)), 1);  med(end + 1) = nan;
if ref.label.classes
    cls = classes_of(v, ref, 0);
    for k = floor(min(ref.value) + 0.5):floor(max(ref.value) + 0.5)
        m = cls == k;
        if k == 0 && ~isempty(ref.label.zero_name)
            items{end + 1} = ref.label.zero_name; %#ok<AGROW>
        else
            items{end + 1} = strrep(ref.label.class_name, '{}', sprintf('%d', k)); %#ok<AGROW>
        end
        px(end + 1) = sum(m(:)); frac(end + 1) = px(end) / total; med(end + 1) = median(res.residual(m)); %#ok<AGROW>
    end
else
    w = nice_step(max(ref.value) - min(ref.value));
    edges = floor(min(ref.value) / w) * w:w:max(ref.value) + w * 0.999;
    for k = 1:numel(edges) - 1
        if k < numel(edges) - 1, m = v >= edges(k) & v < edges(k + 1); else, m = v >= edges(k) & v <= edges(k + 1); end
        items{end + 1} = strtrim(sprintf('%g-%g %s', edges(k), edges(k + 1), unit)); %#ok<AGROW>
        px(end + 1) = sum(m(:)); frac(end + 1) = px(end) / total; med(end + 1) = median(res.residual(m)); %#ok<AGROW>
    end
end
T = table(items(:), px(:), frac(:), med(:), 'VariableNames', {'item', 'pixels', 'fraction', 'median_residual'});
end
