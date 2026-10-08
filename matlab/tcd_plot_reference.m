function fig = tcd_plot_reference(ref, file, visible)
%TCD_PLOT_REFERENCE Reference sheet: the simulated colours of a system and where they repeat.
%
%   tcd_plot_reference(ref)                       % on screen
%   tcd_plot_reference(ref, 'sheet.png', false)   % saved, not shown
%
%   Panels: (1) the simulated colour against the label (thickness, layer number, ...);
%   (2) Delta E from the calibration row, and Delta E to the nearest look-alike, i.e. the
%   closest colour of a structure more than label.gap away - where it drops below ~2 the colour
%   alone cannot tell the two apart; (3) the colour path in the CAM02-UCS a'b' plane, whose
%   loops are the repeating colours; (4) a colour chart of swatches to compare with the eyepiece.
%   tcd_build_reference writes this sheet next to the reference when given 'Out'.
%
%   See also TCD_BUILD_REFERENCE, TCD_SPECTRUM.
if nargin < 2, file = ''; end
if nargin < 3, visible = isempty(file); end
vis = 'off'; if visible, vis = 'on'; end
[x, order] = sort(ref.value);
amb = ref_ambiguity(ref);
steps = diff(x);
typical = median(steps(steps > 0));
if isempty(typical) || isnan(typical), typical = 1; end
unit = ref.label.unit;
axisTitle = ref.label.title;
if ~isempty(unit) && ~contains(lower(axisTitle), regexprep(lower(unit), 's$', ''))
    axisTitle = sprintf('%s (%s)', axisTitle, unit);
end

fig = figure('Visible', vis, 'Position', [60 40 1200 960], 'Color', 'w');

% (1) colour strip, resampled on the label axis (rows need not be evenly spaced)
ax = axes(fig, 'Position', [0.07 0.855 0.88 0.08]);
xs = linspace(x(1), x(end), 1500);
[~, j] = min(abs(x - xs), [], 1);
strip = ref.sRGB(order(j), :);
strip(abs(x(j)' - xs) > 2.5 * typical, :) = 1;   % white where the sweep has no candidates
image(ax, xs, [0 1], reshape(strip, 1, [], 3));
set(ax, 'YTick', [], 'YDir', 'normal');
xlim(ax, [x(1) x(end) + eps]);
ttl = sprintf('%s: simulated colour under %s, NA = %g', ref.name, ref.meta.illuminant, ref.meta.NA);
if isfield(ref.meta, 'description') && ~isempty(ref.meta.description)
    ttl = [ttl '  -  ' ref.meta.description];
end
title(ax, ttl, 'Interpreter', 'none', 'FontWeight', 'normal');

% (2) Delta E curves, broken across sweep gaps
ax = axes(fig, 'Position', [0.07 0.53 0.88 0.26]);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
cut = find(steps > 5 * typical);
[xb, de, am] = deal(x, ref.dE(order), amb(order));
for c = fliplr(cut')
    xb = [xb(1:c); nan; xb(c + 1:end)];
    de = [de(1:c); nan; de(c + 1:end)];
    am = [am(1:c); nan; am(c + 1:end)];
end
patch(ax, [x(1) x(end) x(end) x(1)], [0 0 2 2], [0.77 0.31 0.32], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(ax, xb, de, 'Color', [0.16 0.44 0.59], 'LineWidth', 1.6);
plot(ax, xb, am, 'Color', [0.77 0.31 0.32], 'LineWidth', 1.2);
legend(ax, {'below ~2 \DeltaE: hard to tell apart in a micrograph', '\DeltaE from the calibration row (row 1)', ...
    sprintf('\\DeltaE to the nearest look-alike more than %g %s away', ref.label.gap, unit)}, 'Location', 'best');
xlim(ax, [x(1) x(end) + eps]);
ylim(ax, [0 inf]);
xlabel(ax, axisTitle);
ylabel(ax, '\DeltaE (CAM02-UCS)');

% (3) colour path in a'b'
ax = axes(fig, 'Position', [0.07 0.09 0.36 0.36]);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
[~, ab1, ab2] = deal(0, ref.Jab(order, 2), ref.Jab(order, 3));
for c = fliplr(cut')
    ab1 = [ab1(1:c); nan; ab1(c + 1:end)];
    ab2 = [ab2(1:c); nan; ab2(c + 1:end)];
end
plot(ax, ab1, ab2, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
scatter(ax, ref.Jab(order, 2), ref.Jab(order, 3), 10, x, 'filled');
plot(ax, ref.Jab(1, 2), ref.Jab(1, 3), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = axisTitle;
axis(ax, 'equal');
xlabel(ax, 'a'''); ylabel(ax, 'b''');
title(ax, 'Colour path in CAM02-UCS (loops = repeating colours)', 'FontWeight', 'normal');

% (4) colour chart
ax = axes(fig, 'Position', [0.55 0.09 0.4 0.36]);
rows = chart_rows(ref);
ncol = min(4, numel(rows)); if numel(rows) > 8, ncol = 4; end
nrow = ceil(numel(rows) / ncol);
hold(ax, 'on');
for k = 1:numel(rows)
    r = rows(k);
    cx = mod(k - 1, ncol);  cy = nrow - 1 - floor((k - 1) / ncol);
    rectangle(ax, 'Position', [cx + 0.04, cy + 0.06, 0.92, 0.88], 'FaceColor', ref.sRGB(r, :), 'EdgeColor', 'none');
    if ref.label.classes, txt = class_name(ref, round(ref.value(r))); else, txt = value_text(ref, ref.value(r)); end
    ink = 'w'; if ref.sRGB(r, :) * [0.299; 0.587; 0.114] > 0.5, ink = 'k'; end
    text(ax, cx + 0.5, cy + 0.5, txt, 'HorizontalAlignment', 'center', 'Color', ink, 'Interpreter', 'none');
end
axis(ax, 'equal'); axis(ax, [0 ncol 0 nrow]); axis(ax, 'off');
title(ax, 'Colour chart (sRGB, D65 display)', 'FontWeight', 'normal', 'Visible', 'on');

annotation(fig, 'textbox', [0.07 0.0 0.9 0.035], 'String', ref.meta.stack, 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'FontSize', 7, 'Color', [0.35 0.35 0.35], 'FontName', 'FixedWidth');
if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    if ~visible, close(fig); end
end
end

function rows = chart_rows(ref)
% Every class, or about 13 evenly spaced label values.
if ref.label.classes
    ks = floor(min(ref.value) + 0.5):floor(max(ref.value) + 0.5);
    rows = arrayfun(@(k) ref_find(ref, sprintf('%d', k)), ks);
else
    w = nice_step(max(ref.value) - min(ref.value), 12);
    targets = ceil(min(ref.value) / w) * w:w:max(ref.value) + w * 1e-6;
    [~, rows] = min(abs(ref.value - targets), [], 1);
    rows = unique(rows, 'stable');
end
end

function s = class_name(ref, k)
if k == 0 && ~isempty(ref.label.zero_name)
    s = ref.label.zero_name;
else
    s = strrep(ref.label.class_name, '{}', sprintf('%d', k));
end
end

function s = value_text(ref, v)
s = sprintf('%g', v);
if ~isempty(ref.label.unit) && ~strcmp(ref.label.unit, 'layers'), s = [s ' ' ref.label.unit]; end
end
