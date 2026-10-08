function sys = tcd_load_system(source, opts)
%TCD_LOAD_SYSTEM Read and check a system file: what the sample is made of and what varies.
%
%   sys = tcd_load_system('moo3')                                % bundled: ../systems/moo3.jsonc
%   sys = tcd_load_system('../systems/my_sample.jsonc')
%   sys = tcd_load_system('ps', 'Set', {'oxide_nm', 285; 'bead_nm', 500})   % change constants
%
%   A system file (JSON with // comments, see systems/README.md and systems/template.jsonc)
%   lists the ambient medium, the layers from the light side down, the substrate, and a sweep
%   of structural parameters. Every combination of the sweep is one candidate structure; this
%   function expands them into sys.rows and checks every material and expression, so mistakes
%   show up here rather than halfway through a long calculation.
%
%   The same files are read by the Python code (python/tcd/system.py), with the same rules.
%
%   sys fields:
%     name, description, file       from the file
%     constants                     struct of numbers (after 'Set')
%     rows                          struct of N-by-1 columns: sweep parameters, then the label
%     params                        names of the sweep parameters
%     nRows                         number of candidate structures; row 1 is the calibration row
%     label                         name, unit, title, classes, class_name, zero_name, gap
%     optics                        illuminant, na, na_weighting, angles
%     layers, materials, ambient, substrate, grids, definition   (parsed file contents)
%
%   Options (name-value):
%     Set ({})      constants to change: {'name', value; ...}, {'name', value, ...} or a struct.
%                   A value can be a number or an expression of the other constants.
%
%   See also TCD_BUILD_REFERENCE, TCD_SPECTRUM, TCD_MATERIALS.
arguments
    source
    opts.Set = {}
    opts.Folder (1,:) char = ''     % internal: folder for relative n,k files of a struct definition
end
root = repo_root();
if isstruct(source)
    def = source;
    file = '';
    folder = opts.Folder;
else
    file = find_system_file(char(source), root);
    def = read_jsonc(file);
    folder = fileparts(file);
end
where = file;
if isempty(where), where = 'system'; end
fail = @(fmt, varargin) error('tcd:system', ['%s: ' fmt], where, varargin{:});

check_keys(def, {'name', 'description', 'constants', 'materials', 'ambient', 'layers', 'substrate', ...
    'sweep', 'label', 'optics', 'notes'}, 'top level', where);
for key = {'name', 'layers', 'substrate', 'sweep'}
    if ~isfield(def, key{1}), fail('missing required key ''%s''', key{1}); end
end
sys.name = char(def.name);
sys.description = char(get_or(def, 'description', ''));
sys.file = file;
sys.folder = folder;
sys.where = where;
sys.definition = def;

% constants, in file order, each may use the ones before it
raw = get_or(def, 'constants', struct());
if ~isstruct(raw), fail('''constants'' must be an object of name: value'); end
over = normalise_set(opts.Set);
unknown = setdiff(fieldnames(over), fieldnames(raw));
if ~isempty(unknown)
    fail('cannot set %s: not a constant of this system (constants: %s)', strjoin(unknown, ', '), ...
        strjoin(fieldnames(raw)', ', '));
end
sys.constants = struct();
for k = fieldnames(raw)'
    check_name(k{1}, 'constant', where);
    v = raw.(k{1});
    if isfield(over, k{1}), v = over.(k{1}); end
    try
        sys.constants.(k{1}) = double(expr_eval(v, sys.constants));
    catch err
        fail('constant ''%s'': %s', k{1}, err.message);
    end
end

sys.materials = get_or(def, 'materials', struct());
if ~isstruct(sys.materials), fail('''materials'' must be an object of name: material'); end
for k = fieldnames(sys.materials)'
    check_name(k{1}, 'material', where);
end
sys.ambient = get_or(def, 'ambient', 'Air');
sys.substrate = def.substrate;
sys.layers = normalise_layers(def.layers, 'layers', where);

% optics
o = get_or(def, 'optics', struct());
check_keys(o, {'illuminant', 'na', 'na_weighting', 'angles', 'notes'}, '''optics''', where);
sys.optics = struct('illuminant', char(get_or(o, 'illuminant', 'D65')), ...
    'na', double(expr_eval(get_or(o, 'na', 0), sys.constants)), ...
    'na_weighting', char(get_or(o, 'na_weighting', 'uniform')), 'angles', double(get_or(o, 'angles', 24)));
if sys.optics.na < 0 || sys.optics.na >= 1
    fail('optics.na must be in [0, 1) (objective numerical aperture in air)');
end
if ~any(strcmp(sys.optics.na_weighting, {'uniform', 'gaussian'}))
    fail('optics.na_weighting must be ''uniform'' or ''gaussian''');
end

% sweep -> candidate structures
sys.grids = as_list(def.sweep);
[sys.rows, sys.params] = expand_sweep(sys.grids, sys.constants, where);
sys.nRows = numel(sys.rows.(sys.params{1}));
[sys.label, sys.rows] = parse_label(get_or(def, 'label', struct()), sys.rows, sys.params, sys.constants, where);
sys.vars = sys.constants;
for k = fieldnames(sys.rows)'
    sys.vars.(k{1}) = sys.rows.(k{1});
end
sys.cache = struct('cols', containers.Map(), 'nk', containers.Map());   % handles: shared by copies

% resolve every material once up front so typos surface now
for i = unique([1, floor(sys.nRows / 2) + 1, sys.nRows])
    system_stack(sys, i);
end
end

% ---------------------------------------------------------------------------------------------
function file = find_system_file(name, root)
cands = {name, fullfile(root, 'systems', [name '.jsonc']), fullfile(root, 'systems', [name '.json']), ...
    fullfile(root, 'systems', [lower(name) '.jsonc'])};
for c = cands
    [~, ~, ext] = fileparts(c{1});
    if ~isempty(ext) && isfile(c{1})
        file = c{1};
        return
    end
end
d = dir(fullfile(root, 'systems', '*.jsonc'));
known = regexprep({d.name}, '\.jsonc$', '');
error('tcd:system', 'No system file ''%s'' (bundled systems: %s)', name, strjoin(known, ', '));
end

function v = get_or(s, key, default)
if isstruct(s) && isfield(s, key), v = s.(key); else, v = default; end
end

function check_keys(s, allowed, where_in, where)
if ~isstruct(s), return; end
bad = setdiff(fieldnames(s), allowed);
if ~isempty(bad)
    error('tcd:system', '%s: unknown key(s) %s in %s; allowed: %s', where, strjoin(bad', ', '), ...
        where_in, strjoin(sort(allowed), ', '));
end
end

function check_name(name, what, where)
reserved = {'pi', 'sqrt', 'exp', 'log', 'log10', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'abs', ...
    'floor', 'ceil', 'round', 'min', 'max'};
if isempty(regexp(name, '^[A-Za-z_]\w*$', 'once')) || any(strcmp(name, reserved))
    error('tcd:system', '%s: ''%s'' cannot be a %s name (letters, digits and _, starting with a letter, and not a function name)', ...
        where, name, what);
end
end

function over = normalise_set(s)
if isstruct(s), over = s; return; end
if isempty(s), over = struct(); return; end
if size(s, 2) == 2 && ~isvector(s)
    pairs = s;
elseif mod(numel(s), 2) == 0
    pairs = reshape(s, 2, [])';
else
    error('tcd:system', 'Set must be {''name'', value; ...} or a struct');
end
over = struct();
for i = 1:size(pairs, 1)
    over.(char(pairs{i, 1})) = pairs{i, 2};
end
end

function L = as_list(x)
% JSON arrays of objects decode to struct arrays or cells; return a 1-by-N cell either way.
if iscell(x)
    L = x(:)';
elseif isstruct(x)
    L = num2cell(x(:)');
else
    L = {x};
end
end

function layers = normalise_layers(raw, where_in, where)
layers = as_list(raw);
for j = 1:numel(layers)
    item = layers{j};
    loc = sprintf('%s{%d}', where_in, j);
    if ~isstruct(item)
        error('tcd:system', '%s: %s must be an object with material and thickness', where, loc);
    end
    if isfield(item, 'name'), loc = sprintf('%s (''%s'')', loc, item.name); end
    if isfield(item, 'layers')
        check_keys(item, {'name', 'layers', 'repeat', 'notes'}, loc, where);
        if isempty(item.layers)
            error('tcd:system', '%s: %s: ''layers'' of a group must be a non-empty list', where, loc);
        end
        item.layers = normalise_layers(item.layers, [loc '.layers'], where);
    else
        check_keys(item, {'name', 'material', 'thickness', 'repeat', 'notes'}, loc, where);
        for key = {'material', 'thickness'}
            if ~isfield(item, key{1}), error('tcd:system', '%s: %s needs ''%s''', where, loc, key{1}); end
        end
    end
    layers{j} = item;
end
end

function v = sweep_values(spec, pname, constants, where)
% Swept values of one parameter (column), or [] with derived = true for an expression.
bad = @() error('tcd:system', ['%s: sweep ''%s'': give a number, a list, {"from", "to", "step"}, ' ...
    '{"from", "to", "num"}, {"values": [...]} or an expression string'], where, pname);
if islogical(spec)
    error('tcd:system', '%s: sweep ''%s'': use numbers, not true/false', where, pname);
elseif isnumeric(spec)
    v = double(spec(:));
elseif iscell(spec)
    v = cellfun(@(e) double(expr_eval(e, constants)), spec(:));
elseif isstruct(spec)
    ev = @(key) double(expr_eval(spec.(key), constants));
    if isfield(spec, 'values')
        vals = spec.values;
        if iscell(vals), v = cellfun(@(e) double(expr_eval(e, constants)), vals(:)); else, v = double(vals(:)); end
    elseif all(isfield(spec, {'from', 'to', 'step'}))
        a = ev('from'); b = ev('to'); s = ev('step');
        if s <= 0 || b < a
            error('tcd:system', '%s: sweep ''%s'': need step > 0 and to >= from', where, pname);
        end
        n = floor((b - a) / s + 1e-9) + 1;
        v = round(a + s * (0:n - 1)', 10);
    elseif all(isfield(spec, {'from', 'to', 'num'}))
        v = round(linspace(ev('from'), ev('to'), ev('num'))', 10);
    else
        bad();
    end
else
    bad();
end
end

function [rows, names] = expand_sweep(grids, constants, where)
names = {};
parts = {};
for g = 1:numel(grids)
    blk = grids{g};
    if ~isstruct(blk) || isempty(fieldnames(blk))
        error('tcd:system', '%s: sweep block %d must be an object of parameter: values', where, g);
    end
    if isfield(blk, 'notes'), blk = rmfield(blk, 'notes'); end
    gnames = fieldnames(blk)';
    for k = gnames
        check_name(k{1}, 'sweep parameter', where);
        if isfield(constants, k{1})
            error('tcd:system', '%s: ''%s'' is both a constant and a sweep parameter', where, k{1});
        end
    end
    if g == 1
        names = gnames;
    elseif ~isempty(setxor(gnames, names))
        error('tcd:system', '%s: sweep block %d has parameters %s but block 1 has %s; every block must give every parameter', ...
            where, g, strjoin(sort(gnames), ', '), strjoin(sort(names), ', '));
    end
    derived = cellfun(@(k) ischar(blk.(k)) || isstring(blk.(k)), gnames);
    fixed = gnames(~derived);
    if isempty(fixed)
        error('tcd:system', '%s: sweep block %d has no swept parameter (only expressions)', where, g);
    end
    vals = cellfun(@(k) sweep_values(blk.(k), k, constants, where), fixed, 'UniformOutput', false);
    G = cell(1, numel(fixed));
    [G{end:-1:1}] = ndgrid(vals{end:-1:1});     % first parameter varies slowest, as in Python
    cols = struct();
    for k = 1:numel(fixed)
        cols.(fixed{k}) = G{k}(:);
    end
    n = numel(G{1});
    for k = gnames(derived)                      % derived parameters, in file order
        vars = constants;
        for f = fieldnames(cols)', vars.(f{1}) = cols.(f{1}); end
        try
            cols.(k{1}) = double(expr_eval(blk.(k{1}), vars)) .* ones(n, 1);
        catch err
            error('tcd:system', '%s: sweep parameter ''%s'': %s', where, k{1}, err.message);
        end
    end
    parts{end + 1} = cols; %#ok<AGROW>
end
rows = struct();
for k = names
    rows.(k{1}) = cell2mat(cellfun(@(p) p.(k{1}), parts(:), 'UniformOutput', false));
end
end

function [label, rows] = parse_label(raw, rows, params, constants, where)
check_keys(raw, {'name', 'value', 'unit', 'title', 'classes', 'class_name', 'zero_name', 'gap', 'notes'}, ...
    '''label''', where);
name = get_or(raw, 'name', '');
if isempty(name)
    if numel(params) ~= 1
        error('tcd:system', ['%s: several sweep parameters (%s); say which number the map should report ' ...
            'with label.name (and label.value if it is a combination)'], where, strjoin(params, ', '));
    end
    name = params{1};
end
name = char(name);
check_name(name, 'label', where);
if any(strcmp(name, params))
    if isfield(raw, 'value') && ~strcmp(strtrim(char(string(raw.value))), name)
        error('tcd:system', '%s: label ''%s'' is a sweep parameter; drop label.value or give the label a new name', where, name);
    end
else
    if ~isfield(raw, 'value')
        error('tcd:system', '%s: label ''%s'' is not a sweep parameter, so give label.value, an expression of %s', ...
            where, name, strjoin(params, ', '));
    end
    vars = constants;
    for f = fieldnames(rows)', vars.(f{1}) = rows.(f{1}); end
    try
        rows.(name) = double(expr_eval(raw.value, vars)) .* ones(numel(rows.(params{1})), 1);
    catch err
        error('tcd:system', '%s: label.value: %s', where, err.message);
    end
end
v = rows.(name);
classes = logical(get_or(raw, 'classes', false));
span = max(v) - min(v);
if span == 0, span = 1; end
if classes, gap = 0.5; else, gap = str2double(sprintf('%.3g', span / 15)); end   % 3 significant figures
gap = double(expr_eval(get_or(raw, 'gap', gap), constants));
label = struct('name', name, 'unit', char(get_or(raw, 'unit', '')), 'title', char(get_or(raw, 'title', name)), ...
    'classes', classes, 'class_name', char(get_or(raw, 'class_name', '{}')), ...
    'zero_name', char(get_or(raw, 'zero_name', '')), 'gap', gap);
end
