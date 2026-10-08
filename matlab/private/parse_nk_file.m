function [wl, n, k, meta] = parse_nk_file(file, units)
%PARSE_NK_FILE (wavelength nm, n, k, metadata) from an n,k text file. Same rules as
%   parse_nk_text in python/tcd/materials.py:
%   - columns wavelength, n[, k] separated by commas, tabs, spaces or semicolons;
%   - header lines (not numbers) and '# key: value' metadata lines;
%   - the refractiveindex.info CSV export: an 'wl,n' table followed by an 'wl,k' table;
%   - wavelengths in nm, or um when every value is below 50 (UNITS 'auto'); 'nm' / 'um' force it.
if nargin < 2, units = 'auto'; end
lines = splitlines(string(fileread(file)));
meta = struct();
blocks = struct('header', {}, 'rows', {});
cur = 0;
for i = 1:numel(lines)
    s = strtrim(regexprep(char(lines(i)), '^\xFEFF', ''));
    if isempty(s), continue; end
    if s(1) == '#'
        tok = regexp(s, '^#\s*([A-Za-z][A-Za-z _]*?)\s*:\s*(.*)$', 'tokens', 'once');
        if ~isempty(tok)
            meta.(matlab.lang.makeValidName(lower(strtrim(tok{1})))) = strtrim(tok{2});
        end
        continue
    end
    parts = regexp(s, '[,\s;]+', 'split');
    parts = parts(~cellfun(@isempty, parts));
    vals = str2double(parts);
    if any(isnan(vals))                               % a header line starts a new block
        blocks(end + 1) = struct('header', lower(s), 'rows', {{}}); %#ok<AGROW>
        cur = numel(blocks);
        continue
    end
    if cur == 0
        blocks(end + 1) = struct('header', '', 'rows', {{}}); %#ok<AGROW>
        cur = numel(blocks);
    end
    blocks(cur).rows{end + 1} = vals;
end
blocks = blocks(arrayfun(@(b) ~isempty(b.rows), blocks));
if isempty(blocks)
    error('tcd:nk', '%s: no (wavelength, n[, k]) rows found', file);
end
if numel(blocks) == 1
    d = as_table(blocks(1), file);
    wl = d(:, 1); n = d(:, 2);
    if size(d, 2) > 2, k = d(:, 3); else, k = zeros(size(n)); end
elseif numel(blocks) == 2
    d1 = as_table(blocks(1), file);  d2 = as_table(blocks(2), file);
    isK = arrayfun(@(b) ~isempty(regexp(b.header, '(^|[^a-z])k([^a-z]|$)', 'once')), blocks);
    if isK(1) && ~isK(2), dn = d2; dk = d1; else, dn = d1; dk = d2; end
    wl = dn(:, 1); n = dn(:, 2);
    [kw, o] = sort(dk(:, 1));
    k = interp1(kw, dk(o, 2), wl, 'linear');
    k(wl < kw(1)) = dk(o(1), 2);  k(wl > kw(end)) = dk(o(end), 2);   % held constant outside, as np.interp
else
    error('tcd:nk', '%s: found %d separate tables; expected one (wavelength, n, k) table or an n table followed by a k table', ...
        file, numel(blocks));
end
if strcmp(units, 'um') || (strcmp(units, 'auto') && max(wl) < 50)
    wl = wl * 1000;
elseif ~any(strcmp(units, {'auto', 'nm'}))
    error('tcd:nk', 'units must be ''auto'', ''nm'' or ''um''');
end
[wl, o] = sort(wl);
n = n(o); k = k(o);
if any(diff(wl) == 0)
    error('tcd:nk', '%s: the same wavelength appears twice', file);
end
end

function d = as_table(b, file)
w = min(cellfun(@numel, b.rows));
if w < 2
    error('tcd:nk', '%s: rows need at least two columns (wavelength, n)', file);
end
d = cell2mat(cellfun(@(r) r(1:w), b.rows', 'UniformOutput', false));
end
