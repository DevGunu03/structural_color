function fig = plot_result(res, visible)
%PLOT_RESULT Six-panel figure of a tcd_map_image result.
ref = res.reference;
vis = 'on'; if ~visible, vis = 'off'; end
fig = figure('Visible', vis, 'Position', [50 50 1500 860], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0 1 0.95];                    % leave a strip at the top for the heading
[~, name, ext] = fileparts(res.settings.image);
heading = ref.meta.stack;
if isfield(ref.meta, 'description') && ~isempty(ref.meta.description), heading = ref.meta.description; end
annotation(fig, 'textbox', [0 0.955 1 0.04], 'String', sprintf('%s%s  |  %s: %s', name, ext, ref.name, heading), ...
    'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);

ax = nexttile(tl); image(ax, res.image); axis(ax, 'image', 'off'); hold(ax, 'on');
for k = 1:numel(res.regions)
    r = res.regions(k).roi;
    rectangle(ax, 'Position', [(r(1) - 1) * res.scale + 0.5, (r(2) - 1) * res.scale + 0.5, r(3) * res.scale, r(4) * res.scale], ...
        'EdgeColor', 'w', 'LineWidth', 1.5);
    text(ax, (r(1) - 1) * res.scale, (r(2) - 1) * res.scale - 4, res.regions(k).name, 'Color', 'w', ...
        'BackgroundColor', [0 0 0 0.5], 'FontSize', 8, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end
title(ax, 'Micrograph + calibration regions');

ax = nexttile(tl); image(ax, res.calibrated); axis(ax, 'image', 'off');
title(ax, 'Calibrated to simulated substrate colour');

ax = nexttile(tl); imagesc(ax, res.tcd); axis(ax, 'image', 'off');
v = sort(res.tcd(:)); top = v(max(1, round(0.995 * numel(v))));   % 99.5th percentile, no toolbox
colormap(ax, parula); caxis(ax, [0 max(top, eps)]);
cb = colorbar(ax); cb.Label.String = '\DeltaE from substrate (CAM02-UCS)';
title(ax, 'TCD map (absolute \DeltaE)');

ax = nexttile(tl);
lab = ref.label;
if lab.classes
    ks = floor(min(ref.value) + 0.5):floor(max(ref.value) + 0.5);
    rows = arrayfun(@(k) ref_find(ref, sprintf('%d', k)), ks);
    imagesc(ax, res.classes, [ks(1) - 1.5, ks(end) + 0.5]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; ref.sRGB(rows, :)]);
    names = arrayfun(@(k) strrep(lab.class_name, '{}', sprintf('%d', k)), ks, 'UniformOutput', false);
    if ks(1) == 0 && ~isempty(lab.zero_name), names{1} = lab.zero_name; end
    cb = colorbar(ax); cb.Ticks = ks(1) - 1:ks(end);
    cb.TickLabels = [{'unassigned'}, names];
    title(ax, sprintf('%s (shown in simulated colours)', lab.title), 'Interpreter', 'none');
else
    lo = min(ref.value); hi = max(ref.value);
    v = res.value; v(isnan(v)) = lo - (hi - lo) / 255;
    imagesc(ax, v, [lo - (hi - lo) / 255, hi]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; parula(255)]);
    cb = colorbar(ax); cb.Limits = [lo hi];
    cb.Label.String = lab.title;
    if ~isempty(lab.unit), cb.Label.String = sprintf('%s (%s)', lab.title, lab.unit); end
    cb.Label.Interpreter = 'none';
    title(ax, sprintf('%s (nearest reference, grey = unassigned)', lab.title), 'Interpreter', 'none');
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
