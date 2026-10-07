function fig = plot_result(res, visible)
%PLOT_RESULT Six-panel figure of a tcd_map_image result.
ref = res.reference;
vis = 'on'; if ~visible, vis = 'off'; end
fig = figure('Visible', vis, 'Position', [50 50 1500 860], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0 1 0.95];                    % leave a strip at the top for the heading
[~, name, ext] = fileparts(res.settings.image);
annotation(fig, 'textbox', [0 0.955 1 0.04], 'String', sprintf('%s%s  |  %s', name, ext, ref.meta.stack), ...
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
if strcmp(ref.system, 'PS')
    n = floor(max(ref.value));
    dense = 1;
    for k = 1:n, dense(end + 1) = find(ref.params.layers == k & abs(ref.params.packing - 1) < 1e-9, 1); end %#ok<AGROW>
    imagesc(ax, res.layers, [-1.5, n + 0.5]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; ref.sRGB(dense, :)]);
    cb = colorbar(ax); cb.Ticks = -1:n;
    cb.TickLabels = [{'unassigned', 'substrate'}, arrayfun(@(k) sprintf('%dL', k), 1:n, 'UniformOutput', false)];
    title(ax, 'Layer number (shown in simulated colours)');
else
    v = res.value; v(isnan(v)) = -1;
    imagesc(ax, v, [-max(ref.value) / 255, max(ref.value)]); axis(ax, 'image', 'off');
    colormap(ax, [0.82 0.82 0.82; parula(255)]);
    cb = colorbar(ax); cb.Label.String = 'MoO_3 thickness (nm)'; cb.Limits = [0 max(ref.value)];
    title(ax, 'Thickness (nearest reference, grey = unassigned)');
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
