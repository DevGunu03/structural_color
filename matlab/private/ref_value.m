function v = ref_value(ref, spec)
%REF_VALUE Label value of a text such as '230nm', '230 nm', '2L', '2' or 'substrate'.
%   Same as Reference.parse_value in Python.
low = lower(strrep(strtrim(char(spec)), ' ', ''));
zero = 'substrate';
if isfield(ref.label, 'zero_name') && ~isempty(ref.label.zero_name), zero = lower(ref.label.zero_name); end
if any(strcmp(low, {'substrate', 'bare', zero}))
    v = ref.value(1);
    return
end
for suffix = {lower(strrep(ref.label.unit, ' ', '')), 'nm', 'l'}
    if ~isempty(suffix{1}) && endsWith(low, suffix{1})
        low = low(1:end - numel(suffix{1}));
        break
    end
end
v = str2double(low);
if isnan(v)
    error('tcd:value', 'Cannot read ''%s'': use ''substrate'' or a %s value such as ''250nm'' or ''2L''', spec, ref.label.title);
end
end
