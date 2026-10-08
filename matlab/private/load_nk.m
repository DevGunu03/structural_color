function n = load_nk(material, wl)
%LOAD_NK Complex refractive index n + ik of MATERIAL on wavelengths WL (nm).
%   MATERIAL is a library name, a path to an n,k file, or a number. The library is
%   data/nk_library.csv (columns <name>_n, <name>_k) plus one file per added material in
%   data/materials/ (tcd_add_material); case, '-', '_' and spaces in names are ignored.
%   Linear interpolation and extrapolation, as in TransferMatrix.m.
%   load_nk('list') returns the library as a struct array (name, origin, source, wl, n, k).
persistent lib
if ischar(material) && strcmp(material, 'clear-cache')
    lib = [];
    return
end
if isempty(lib), lib = read_library(); end
if ischar(material) && strcmp(material, 'list')
    v = values(lib);
    n = [v{:}];
    [~, o] = sort(lower({n.name}));
    n = n(o);
    return
end
if isnumeric(material)
    n = complex(material) * ones(size(wl));
    return
end
[~, ~, ext] = fileparts(material);
if any(strcmpi(ext, {'.csv', '.txt', '.dat', '.nk', '.tsv'})) && isfile(material)
    [w, nn, kk] = parse_nk_file(material);
else
    key = normkey(material);
    if ~isKey(lib, key)
        v = values(lib);
        names = cellfun(@(m) m.name, v, 'UniformOutput', false);
        error('tcd:material', 'Material "%s" not in the n,k library. Available: %s', material, strjoin(sort(names), ', '));
    end
    m = lib(key);
    w = m.wl; nn = m.n; kk = m.k;
end
n = interp1(w, nn, wl, 'linear', 'extrap') + 1i * interp1(w, kk, wl, 'linear', 'extrap');
end

function lib = read_library()
root = repo_root();
t = readtable(fullfile(root, 'data', 'nk_library.csv'), 'VariableNamingRule', 'preserve');
names = t.Properties.VariableNames;
w = t{:, 1};
lib = containers.Map();
for c = 2:numel(names)
    nm = names{c};
    if endsWith(nm, '_n')
        base = nm(1:end-2);
        kc = find(strcmp(names, [base '_k']), 1);
        if ~isempty(kc)
            lib(normkey(base)) = struct('name', base, 'origin', 'data/nk_library.csv', ...
                'source', 'nk_library.csv (Index_of_Refraction_library.xls)', 'wl', w, 'n', t{:, c}, 'k', t{:, kc});
        end
    end
end
d = dir(fullfile(root, 'data', 'materials', '*.*'));
for i = 1:numel(d)
    [~, stem, ext] = fileparts(d(i).name);
    if d(i).isdir || ~any(strcmpi(ext, {'.csv', '.txt', '.dat', '.nk', '.tsv'})) || isKey(lib, normkey(stem))
        continue                                 % a built-in name always wins
    end
    try
        [ww, nn, kk, meta] = parse_nk_file(fullfile(d(i).folder, d(i).name));
    catch
        continue
    end
    name = stem;  if isfield(meta, 'name'), name = meta.name; end
    src = '';     if isfield(meta, 'source'), src = meta.source; end
    lib(normkey(stem)) = struct('name', name, 'origin', ['data/materials/' d(i).name], 'source', src, ...
        'wl', ww, 'n', nn, 'k', kk);
end
end

function k = normkey(s)
k = lower(regexprep(char(s), '[^A-Za-z0-9]', ''));
end
