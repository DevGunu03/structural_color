function [N, d, names] = system_stack(sys, i, wl)
%SYSTEM_STACK Complex indices and thicknesses of candidate structure I of a loaded system.
%   [N, d, names] = system_stack(sys, i)        on 360:830 nm
%   N: M-by-W (ambient first, substrate last), d: 1-by-M nm, names: 1-by-M cellstr.
%   Layers of zero thickness (or repeat 0) are left out: they change nothing (I_ab I_bc = I_ac).
%   Mirrors System.stack in python/tcd/system.py.
if nargin < 3, wl = 360:830; end
items = expand(sys, sys.layers, i, 'layers');
M = size(items, 1) + 2;
N = complex(zeros(M, numel(wl)));
N(1, :) = nk_of(sys, concrete(sys, sys.ambient, i, 'ambient', {}), wl);
for j = 1:size(items, 1)
    N(j + 1, :) = nk_of(sys, concrete(sys, items{j, 1}.material, i, sprintf('layer ''%s''', items{j, 3}), {}), wl);
end
N(M, :) = nk_of(sys, concrete(sys, sys.substrate, i, 'substrate', {}), wl);
d = [0, cell2mat(items(:, 2))', 0];
names = [{'ambient'}, items(:, 3)', {'substrate'}];
end

% ---------------------------------------------------------------------------------------------
function out = expand(sys, layers, i, where)
% Rows of {layer struct, thickness nm, name} after applying repeat and dropping zero thickness.
out = cell(0, 3);
for j = 1:numel(layers)
    item = layers{j};
    loc = sprintf('%s{%d}', where, j);
    name = loc;
    if isfield(item, 'name'), name = char(item.name); loc = sprintf('%s (''%s'')', loc, name); end
    rep = 1;
    if isfield(item, 'repeat'), rep = col(sys, item.repeat, [loc ' repeat'], i); end
    if rep < -1e-9 || abs(rep - round(rep)) > 1e-6
        error('tcd:system', '%s: %s: repeat = %g for %s; it must be a whole number >= 0', ...
            sys.where, loc, rep, row_text(sys, i));
    end
    rep = round(rep);
    if isfield(item, 'layers')
        for r = 1:rep
            out = [out; expand(sys, item.layers, i, sprintf('%s{%d}.layers', where, j))]; %#ok<AGROW>
        end
        continue
    end
    t = col(sys, item.thickness, [loc ' thickness'], i);
    if ~isfinite(t) || t < 0
        error('tcd:system', '%s: %s: thickness = %g nm for %s', sys.where, loc, t, row_text(sys, i));
    end
    if t > 0
        out = [out; repmat({item, t, name}, rep, 1)]; %#ok<AGROW>
    end
end
end

function v = col(sys, expr, what, i)
% Value of a number/expression for row i; whole columns are computed once and cached.
if isnumeric(expr) || islogical(expr)
    v = double(expr(1));
    return
end
key = char(expr);
if ~isKey(sys.cache.cols, key)
    try
        c = double(expr_eval(key, sys.vars));
    catch err
        error('tcd:system', '%s: %s: %s', sys.where, what, err.message);
    end
    sys.cache.cols(key) = c .* ones(sys.nRows, 1);
end
c = sys.cache.cols(key);
v = c(i);
end

function s = row_text(sys, i)
f = fieldnames(sys.rows)';
s = strjoin(cellfun(@(k) sprintf('%s=%g', k, sys.rows.(k)(i)), f, 'UniformOutput', false), ', ');
end

function c = concrete(sys, spec, i, where, seen)
% The material of row i with every expression replaced by its number; c.key identifies it.
if isnumeric(spec) && isscalar(spec)
    c = struct('kind', 'n', 'n', double(spec), 'k', 0);
    c.key = sprintf('n:%.17g:%.17g', c.n, c.k);
    return
end
if isstring(spec), spec = char(spec); end
if ischar(spec)
    if isfield(sys.materials, spec)
        if any(strcmp(seen, spec))
            error('tcd:system', '%s: material ''%s'' refers to itself', sys.where, spec);
        end
        c = concrete(sys, sys.materials.(spec), i, sprintf('material ''%s''', spec), [seen, {spec}]);
        return
    end
    [~, ~, ext] = fileparts(spec);
    if any(strcmpi(ext, {'.csv', '.txt', '.dat', '.nk', '.tsv'}))
        c = file_concrete(sys, spec, where);
        return
    end
    try
        load_nk(spec, 500);
    catch err
        error('tcd:system', '%s: %s: ''%s'' is not defined under ''materials'' and %s', sys.where, where, spec, ...
            regexprep(err.message, '^Material "[^"]*" ', ''));
    end
    c = struct('kind', 'library', 'name', spec, 'key', ['lib:' spec]);
    return
end
if ~isstruct(spec)
    error('tcd:system', '%s: %s: a material is a name, a number or an object', sys.where, where);
end
if isfield(spec, 'notes'), spec = rmfield(spec, 'notes'); end
kinds = {'library', {'library'}; 'file', {'file'}; 'n', {'n', 'k'}; 'cauchy', {'cauchy', 'k'}; ...
    'sellmeier', {'sellmeier', 'k'}; 'ema', {'ema', 'host', 'inclusion', 'fraction'}};
kind = find(cellfun(@(k) isfield(spec, k), kinds(:, 1)), 1);
if isempty(kind)
    error('tcd:system', '%s: %s: material object needs one of %s', sys.where, where, strjoin(sort(kinds(:, 1))', ', '));
end
bad = setdiff(fieldnames(spec), kinds{kind, 2});
if ~isempty(bad)
    error('tcd:system', '%s: unknown key(s) %s in %s; allowed: %s', sys.where, strjoin(bad', ', '), where, ...
        strjoin(sort(kinds{kind, 2}), ', '));
end
num = @(key, default) getnum(sys, spec, key, default, where, i);
switch kinds{kind, 1}
    case 'library'
        c = concrete(sys, char(spec.library), i, where, seen);
    case 'file'
        c = file_concrete(sys, char(spec.file), where);
    case 'n'
        c = struct('kind', 'n', 'n', num('n', []), 'k', num('k', 0));
        c.key = sprintf('n:%.17g:%.17g', c.n, c.k);
    case 'cauchy'
        cc = spec.cauchy;
        if isstruct(cc)
            cc = cellfun(@(k) cc.(k), intersect({'A', 'B', 'C'}, fieldnames(cc), 'stable'), 'UniformOutput', false);
        elseif isnumeric(cc)
            cc = num2cell(cc(:)');
        end
        if ~iscell(cc) || isempty(cc) || numel(cc) > 3
            error('tcd:system', '%s: %s: cauchy is [A, B, C] or {"A":..,"B":..,"C":..}', sys.where, where);
        end
        coef = cellfun(@(e) col(sys, e, [where '.cauchy'], i), cc);
        c = struct('kind', 'cauchy', 'coef', coef, 'k', num('k', 0));
        c.key = ['cauchy:' sprintf('%.17g:', coef, c.k)];
    case 'sellmeier'
        s = spec.sellmeier;
        if ~isstruct(s) || ~isempty(setxor(fieldnames(s), {'B', 'C'}))
            error('tcd:system', '%s: %s: sellmeier is {"B": [...], "C": [...]} (C in um^2)', sys.where, where);
        end
        B = cellfun(@(e) col(sys, e, [where '.sellmeier.B'], i), to_cell(s.B));
        C = cellfun(@(e) col(sys, e, [where '.sellmeier.C'], i), to_cell(s.C));
        if numel(B) ~= numel(C)
            error('tcd:system', '%s: %s: Sellmeier B and C need the same number of terms', sys.where, where);
        end
        c = struct('kind', 'sellmeier', 'B', B, 'C', C, 'k', num('k', 0));
        c.key = ['sellmeier:' sprintf('%.17g:', B, C, c.k)];
    case 'ema'
        rule = lower(char(spec.ema));
        if ~any(strcmp(rule, {'bruggeman', 'linear', 'maxwell-garnett'}))
            error('tcd:system', '%s: %s: ema must be one of bruggeman, linear, maxwell-garnett', sys.where, where);
        end
        f = num('fraction', []);
        if f < -1e-12 || f > 1 + 1e-12
            error('tcd:system', '%s: %s: fraction = %g for row %d (%s); it must stay within [0, 1]', ...
                sys.where, where, f, i, row_text(sys, i));
        end
        c = struct('kind', 'ema', 'rule', rule, ...
            'host', concrete(sys, spec.host, i, [where '.host'], seen), ...
            'inclusion', concrete(sys, spec.inclusion, i, [where '.inclusion'], seen), 'f', min(max(f, 0), 1));
        c.key = sprintf('ema:%s:(%s):(%s):%.17g', rule, c.host.key, c.inclusion.key, c.f);
end
end

function c = to_cell(x)
% numeric vector or cell of numbers/expressions -> 1-by-N cell
if iscell(x), c = x(:)'; else, c = num2cell(x(:)'); end
end

function v = getnum(sys, spec, key, default, where, i)
if ~isfield(spec, key)
    if isempty(default)
        error('tcd:system', '%s: %s: missing ''%s''', sys.where, where, key);
    end
    v = default;
    return
end
v = col(sys, spec.(key), [where '.' key], i);
end

function c = file_concrete(sys, name, where)
root = repo_root();
if ~isempty(regexp(name, '^([A-Za-z]:[\\/]|[\\/])', 'once'))   % absolute path
    cands = {name};
else
    cands = {fullfile(root, name), fullfile(root, 'data', 'materials', name)};
    if ~isempty(sys.folder), cands = [{fullfile(sys.folder, name)}, cands]; end
end
hit = find(cellfun(@isfile, cands), 1);
if isempty(hit)
    error('tcd:system', '%s: %s: n,k file ''%s'' not found (looked next to the system file, in the repository root and in data/materials/)', ...
        sys.where, where, name);
end
c = struct('kind', 'file', 'name', cands{hit}, 'key', ['file:' cands{hit}]);
end

function n = nk_of(sys, c, wl)
% n + ik of a concrete material, cached by its key.
key = sprintf('%s@%d', c.key, numel(wl));
if isKey(sys.cache.nk, key)
    n = sys.cache.nk(key);
    return
end
lam = wl / 1000;                                 % um, for the dispersion formulas
switch c.kind
    case 'n'
        n = complex(c.n, c.k) * ones(size(wl));
    case {'library', 'file'}
        n = load_nk(c.name, wl);
    case 'cauchy'
        A = [c.coef, 0, 0];
        n = A(1) + A(2) ./ lam .^ 2 + A(3) ./ lam .^ 4 + 1i * c.k;
    case 'sellmeier'
        e = ones(size(wl));
        for t = 1:numel(c.B)
            e = e + c.B(t) * lam .^ 2 ./ (lam .^ 2 - c.C(t));
        end
        n = real(sqrt(complex(e))) + 1i * c.k;
    case 'ema'
        n = ema_mix(c.rule, nk_of(sys, c.host, wl), nk_of(sys, c.inclusion, wl), c.f);
end
sys.cache.nk(key) = n;
end
