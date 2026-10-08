function row = ref_find(ref, spec)
%REF_FIND Row of a reference for an anchor label (same rules as Reference.find in Python).
%   'substrate' / 'bare' (or the label's zero_name) -> row 1
%   a label value: '2', '2L', '250nm', '250 nm'  -> the row whose label is nearest
%   exact parameters: 'layers=2,packing=0.5'       -> the first row with those values
s = strtrim(char(spec));
low = lower(strrep(s, ' ', ''));
zero = 'substrate';
if isfield(ref.label, 'zero_name') && ~isempty(ref.label.zero_name), zero = lower(ref.label.zero_name); end
if any(strcmp(low, {'substrate', 'bare', zero}))
    row = 1;
    return
end
if contains(s, '=')
    mask = true(numel(ref.value), 1);
    for part = strsplit(s, ',')
        kv = strtrim(strsplit(part{1}, '='));
        if ~isfield(ref.params, kv{1})
            error('tcd:anchor', '''%s'' is not a column of this reference (%s)', kv{1}, strjoin(fieldnames(ref.params)', ', '));
        end
        mask = mask & abs(ref.params.(kv{1}) - str2double(kv{2})) <= 1e-6;
    end
    row = find(mask, 1);
    if isempty(row), error('tcd:anchor', 'No reference row with %s', s); end
    return
end
num = low;
for suffix = {lower(strrep(ref.label.unit, ' ', '')), 'nm', 'l'}
    if ~isempty(suffix{1}) && endsWith(num, suffix{1})
        num = num(1:end - numel(suffix{1}));
        break
    end
end
target = str2double(num);
if isnan(target)
    error('tcd:anchor', 'Cannot read anchor ''%s'': use ''substrate'', a %s value or ''param=value,...''', s, ref.label.title);
end
[~, row] = min(abs(ref.value - target));
end
