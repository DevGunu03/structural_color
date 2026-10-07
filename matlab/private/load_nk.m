function n = load_nk(material, wl)
%LOAD_NK Complex refractive index n + ik of MATERIAL on wavelengths WL (nm).
%   MATERIAL is a column name of data/nk_library.csv without the _n/_k suffix
%   (e.g. 'SiO2-Franta', 'aMoO3', 'PS-beads'; case and -/_ are ignored), a path to a
%   text file of (wavelength, n[, k]) - wavelength in um if all values are < 50 - or a
%   number. Linear interpolation and extrapolation, as in TransferMatrix.m.
persistent lib
if isnumeric(material)
    n = complex(material) * ones(size(wl));
    return
end
if isfile(material)
    d = readmatrix(material, 'FileType', 'text');
    d = d(all(isfinite(d(:, 1:2)), 2), :);
    w = d(:, 1);
    if max(w) < 50, w = w * 1000; end
    nn = d(:, 2);
    if size(d, 2) > 2, kk = d(:, 3); else, kk = zeros(size(nn)); end
    [w, o] = sort(w);
    nn = nn(o); kk = kk(o);
else
    if isempty(lib), lib = read_library(); end
    key = normkey(material);
    if ~isKey(lib, key)
        error('tcd:material', 'Material "%s" not in data/nk_library.csv. Available: %s', ...
            material, strjoin(sort(keys(lib)), ', '));
    end
    d = lib(key);
    w = d(:, 1); nn = d(:, 2); kk = d(:, 3);
end
n = interp1(w, nn, wl, 'linear', 'extrap') + 1i * interp1(w, kk, wl, 'linear', 'extrap');
end

function lib = read_library()
f = fullfile(repo_root(), 'data', 'nk_library.csv');
t = readtable(f, 'VariableNamingRule', 'preserve');
names = t.Properties.VariableNames;
w = t{:, 1};
lib = containers.Map();
for c = 2:numel(names)
    nm = names{c};
    if endsWith(nm, '_n')
        base = nm(1:end-2);
        kc = find(strcmp(names, [base '_k']), 1);
        if ~isempty(kc)
            lib(normkey(base)) = [w, t{:, c}, t{:, kc}];
        end
    end
end
end

function k = normkey(s)
k = lower(regexprep(char(s), '[^A-Za-z0-9]', ''));
end
