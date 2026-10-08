function S = tcd_spectrum(system, at, opts)
%TCD_SPECTRUM Reflectance spectrum and colour of chosen structures of a system file.
%
%   tcd_spectrum('moo3', struct('thickness_nm', [100 200 300]))          % plot three flakes
%   S = tcd_spectrum('ps', struct('layers', [1 2], 'packing', [1 0.5]))  % values, no plot
%   tcd_spectrum('../systems/my_sample.jsonc')                           % first, middle, last row
%
%   Use it to check a new system file before building the whole reference: are the layers in
%   the right order, are the thicknesses what you meant, do the colours look plausible?
%   AT is a struct with one field per sweep parameter (derived parameters are computed); each
%   field holds one value per structure (or one value used for all of them).
%
%   Options (name-value):
%     Set ({})          change constants: {'oxide_nm', 285}
%     Illuminant ('')   override optics.illuminant
%     NA ([])           override optics.na
%     Plot ([])         draw the figure (default: when no output is requested)
%     Out ('')          save the figure to this file
%
%   S fields: wl (1-by-471), R (rows x 471), XYZ, Jab, sRGB, dE (from the first structure),
%   rows (struct of parameters), names (cellstr), stacks (cell of {layer names; thicknesses}).
%
%   See also TCD_LOAD_SYSTEM, TCD_BUILD_REFERENCE, TCD_TMM.
arguments
    system
    at = struct()
    opts.Set = {}
    opts.Illuminant (1,:) char = ''
    opts.NA double = []
    opts.Plot = []
    opts.Out (1,:) char = ''
end
if isstruct(system) && isfield(system, 'rows')
    sys = system;
else
    sys = tcd_load_system(system, 'Set', opts.Set);
end
g0 = sys.grids{1};
if isfield(g0, 'notes'), g0 = rmfield(g0, 'notes'); end
names = fieldnames(g0)';
derived = cellfun(@(k) ischar(g0.(k)) || isstring(g0.(k)), names);
needed = names(~derived);
if isempty(fieldnames(at))                        % default: first, middle and last candidate
    pick = unique([1, floor(sys.nRows / 2) + 1, sys.nRows]);
    at = struct();
    for k = needed, at.(k{1}) = sys.rows.(k{1})(pick); end
end
if ~isempty(setxor(fieldnames(at), needed))
    error('tcd:spectrum', 'Give exactly the parameters %s (got %s)', strjoin(needed, ', '), strjoin(fieldnames(at)', ', '));
end
n = max(structfun(@numel, at));
grids = cell(1, n);
for i = 1:n
    g = struct();
    for k = names
        if any(strcmp(k{1}, needed))
            v = at.(k{1});
            g.(k{1}) = v(min(i, numel(v)));
        else
            g.(k{1}) = g0.(k{1});
        end
    end
    grids{i} = g;
end
def = sys.definition;
def.constants = sys.constants;
def.sweep = grids;
if ~isfield(def, 'label'), def.label = struct(); end
if ~isfield(def.label, 'gap'), def.label.gap = sys.label.gap; end
sub = tcd_load_system(def, 'Folder', sys.folder);

illuminant = opts.Illuminant;
if isempty(illuminant), illuminant = sys.optics.illuminant; end
na = opts.NA;
if isempty(na), na = sys.optics.na; end
wl = 360:830;
S.wl = wl;
S.R = zeros(sub.nRows, numel(wl));
S.stacks = cell(sub.nRows, 1);
S.names = cell(sub.nRows, 1);
for i = 1:sub.nRows
    [N, d, lnames] = system_stack(sub, i, wl);
    S.R(i, :) = tcd_tmm(N, d, wl, 'NA', na, 'Angles', sys.optics.angles, 'Weighting', sys.optics.na_weighting);
    S.stacks{i} = {lnames; d};
    f = fieldnames(sub.rows)';
    S.names{i} = strjoin(cellfun(@(k) sprintf('%s=%g', k, sub.rows.(k)(i)), f, 'UniformOutput', false), ', ');
end
c = colour_const();
S.XYZ = refl_to_xyz(S.R, illuminant);
S.Jab = xyz_to_cam02ucs(S.XYZ);
S.sRGB = min(max(srgb_codec(S.XYZ / 100 * c.XYZ2RGB', 'encode'), 0), 1);
S.dE = vecnorm(S.Jab - S.Jab(1, :), 2, 2);
S.rows = sub.rows;

doPlot = opts.Plot;
if isempty(doPlot), doPlot = nargout == 0 || ~isempty(opts.Out); end
if doPlot
    vis = 'on'; if nargout > 0 && ~isempty(opts.Out), vis = 'off'; end
    fig = figure('Visible', vis, 'Position', [80 80 1100 440], 'Color', 'w');
    ax = axes(fig, 'Position', [0.07 0.13 0.62 0.75]);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    for i = 1:sub.nRows
        plot(ax, wl, S.R(i, :), 'LineWidth', 1.6, 'Color', min(S.sRGB(i, :) * 0.85, 1));
    end
    xlabel(ax, 'wavelength (nm)'); ylabel(ax, 'reflectance'); xlim(ax, [wl(1) wl(end)]); ylim(ax, [0 inf]);
    legend(ax, S.names, 'Interpreter', 'none', 'Location', 'best');
    title(ax, sprintf('%s: %s', sys.name, system_describe(sub)), 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
    ax = axes(fig, 'Position', [0.74 0.13 0.22 0.75]);
    for i = 1:sub.nRows
        y = sub.nRows - i;
        rectangle(ax, 'Position', [0 y + 0.08 1 0.84], 'FaceColor', S.sRGB(i, :), 'EdgeColor', 'none');
        ink = 'w'; if S.sRGB(i, :) * [0.299; 0.587; 0.114] > 0.5, ink = 'k'; end
        text(ax, 0.5, y + 0.5, S.names{i}, 'HorizontalAlignment', 'center', 'Color', ink, 'Interpreter', 'none', 'FontSize', 8);
    end
    axis(ax, [0 1 0 sub.nRows]); axis(ax, 'off');
    if ~isempty(opts.Out)
        exportgraphics(fig, opts.Out, 'Resolution', 150);
        if strcmp(vis, 'off'), close(fig); end
    end
end
if nargout == 0
    for i = 1:sub.nRows
        fprintf('%s\n    stack: %s\n    XYZ = [%s], J''a''b'' = [%s], sRGB = [%s], dE from first = %.2f\n', S.names{i}, ...
            strjoin(cellfun(@(nm, t) strtrim(sprintf('%s %s', nm, num2str_nm(t))), S.stacks{i}{1}, num2cell(S.stacks{i}{2}), ...
            'UniformOutput', false), ' | '), num2str(S.XYZ(i, :), '%.3f '), num2str(S.Jab(i, :), '%.2f '), ...
            num2str(S.sRGB(i, :), '%.3f '), S.dE(i));
    end
    clear S
end
end

function s = num2str_nm(t)
if t == 0, s = ''; else, s = sprintf('%g nm', t); end
end
