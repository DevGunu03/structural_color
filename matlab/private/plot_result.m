function fig = plot_result(res, visible)
%PLOT_RESULT Eight-panel figure of a tcd_map_image result (same panels as the Python figure).
ref = res.reference;
vis = 'on'; if ~visible, vis = 'off'; end
fig = figure('Visible', vis, 'Position', [30 50 1900 860], 'Color', 'w');
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0 1 0.95];                    % leave a strip at the top for the heading
[~, name, ext] = fileparts(res.settings.image);
heading = ref.meta.stack;
if isfield(ref.meta, 'description') && ~isempty(ref.meta.description), heading = ref.meta.description; end
annotation(fig, 'textbox', [0 0.955 1 0.04], 'String', sprintf('%s%s  |  %s: %s', name, ext, ref.name, heading), ...
    'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);

% 1. micrograph with calibration regions (white) and checks (green = agrees, red = disagrees)
ax = nexttile(tl); image(ax, res.image); axis(ax, 'image', 'off'); hold(ax, 'on');
for k = 1:numel(res.regions)
    draw_box(ax, res.regions(k).roi, res.scale, [1 1 1], res.regions(k).name, [0 0 0]);
end
for k = 1:height(res.checks)
    c = res.checks(k, :);
    ok = c.within_tolerance >= 0.5;
    col = [0.8 0 0]; if ok, col = [0 0.8 0]; end
    draw_box(ax, c.roi, res.scale, col, sprintf('check %s: %.4g', c.check{1}, c.median_value), col * 0.6);
end
title(ax, 'Micrograph + calibration regions');

ax = nexttile(tl); image(ax, res.calibrated); axis(ax, 'image', 'off');
title(ax, 'Calibrated to simulated substrate colour');

ax = nexttile(tl); imagesc(ax, res.tcd); axis(ax, 'image', 'off');
v = sort(res.tcd(:)); top = v(max(1, round(0.995 * numel(v))));   % 99.5th percentile, no toolbox
colormap(ax, parula); caxis(ax, [0 max(top, eps)]);
cb = colorbar(ax); cb.Label.String = '\DeltaE from substrate (CAM02-UCS)';
title(ax, 'TCD map (absolute \DeltaE)');

lab = ref.label;
ax = nexttile(tl);
if lab.classes
    label_map(ax, res.classes, ref, sprintf('%s (simulated colours)', lab.title));
else
    label_map(ax, res.value, ref, sprintf('%s (grey = unassigned)', lab.title));
end

% 5. reliability
ax = nexttile(tl);
rel = res.reliability; r = rel; r(isnan(r)) = -1 / 254;
imagesc(ax, r, [-1 / 254, 1]); axis(ax, 'image', 'off');
colormap(ax, [0.82 0.82 0.82; red_yellow_green(254)]);
cb = colorbar(ax); cb.Limits = [0 1]; cb.Label.String = 'reliability (fit \times uniqueness \times consistency)';
title(ax, sprintf('Reliability (colour error %.1f \\DeltaE)', res.sigma.total));

% 6. value only where reliable
ax = nexttile(tl);
shown = res.value;
shown(~(rel >= res.settings.MinReliability)) = nan;
if lab.classes
    cls = floor(shown + 0.5); cls(isnan(cls)) = -1;
    label_map(ax, cls, ref, sprintf('Only where reliability \\geq %g', res.settings.MinReliability));
else
    label_map(ax, shown, ref, sprintf('Only where reliability \\geq %g', res.settings.MinReliability));
end

ax = nexttile(tl); imagesc(ax, res.residual, [0 max(1.5 * res.settings.MaxResidual, 1)]); axis(ax, 'image', 'off');
colormap(ax, hot); cb = colorbar(ax); cb.Label.String = 'distance to nearest reference (\DeltaE)';
title(ax, sprintf('Match residual (unassigned above %g)', res.settings.MaxResidual));

ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
J = reshape(res.Jab, [], 3);
step = max(1, floor(size(J, 1) / 40000));
scatter(ax, J(1:step:end, 2), J(1:step:end, 3), 2, [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha', 0.15);
scatter(ax, ref.Jab(:, 2), ref.Jab(:, 3), 8, ref.value, 'filled');
colormap(ax, parula);
plot(ax, ref.Jab(1, 2), ref.Jab(1, 3), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
axis(ax, 'equal'); xlabel(ax, 'a'''); ylabel(ax, 'b''');
legend(ax, {'image pixels', 'simulated locus', 'substrate'}, 'Location', 'southeast');
title(ax, 'Image colours vs simulated locus');
end

function draw_box(ax, r, scale, col, txt, bg)
rectangle(ax, 'Position', [(r(1) - 1) * scale + 0.5, (r(2) - 1) * scale + 0.5, r(3) * scale, r(4) * scale], ...
    'EdgeColor', col, 'LineWidth', 1.5);
text(ax, (r(1) - 1) * scale, (r(2) - 1) * scale - 4, txt, 'Color', 'w', 'BackgroundColor', [bg 0.6], ...
    'FontSize', 8, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end

function label_map(ax, values, ref, ttl)
% Classes in simulated colours, or continuous values on parula; grey = not shown.
lab = ref.label;
if lab.classes
    ks = floor(min(ref.value) + 0.5):floor(max(ref.value) + 0.5);
    rows = arrayfun(@(k) ref_find(ref, sprintf('%d', k)), ks);
    v = values; v(isnan(v)) = -1;
    imagesc(ax, v, [ks(1) - 1.5, ks(end) + 0.5]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; ref.sRGB(rows, :)]);
    names = arrayfun(@(k) strrep(lab.class_name, '{}', sprintf('%d', k)), ks, 'UniformOutput', false);
    if ks(1) == 0 && ~isempty(lab.zero_name), names{1} = lab.zero_name; end
    cb = colorbar(ax); cb.Ticks = ks(1) - 1:ks(end);
    cb.TickLabels = [{'n/a'}, names];
else
    lo = min(ref.value); hi = max(ref.value);
    v = values; v(isnan(v)) = lo - (hi - lo) / 255;
    imagesc(ax, v, [lo - (hi - lo) / 255, hi]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; parula(255)]);
    cb = colorbar(ax); cb.Limits = [lo hi];
    cb.Label.String = lab.title;
    if ~isempty(lab.unit), cb.Label.String = sprintf('%s (%s)', lab.title, lab.unit); end
    cb.Label.Interpreter = 'none';
end
title(ax, ttl, 'Interpreter', 'tex');
end

function c = red_yellow_green(n)
% Red (0) -> yellow (0.5) -> green (1), as matplotlib's RdYlGn ends.
t = linspace(0, 1, n)';
c = [min(1, 2 * (1 - t)) * 0.84 + 0.1, min(1, 2 * t) * 0.65 + 0.1, 0.15 * ones(n, 1)];
end
