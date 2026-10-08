function tcd_save_reference(ref, file, opts)
%TCD_SAVE_REFERENCE Write a reference as CSV + JSON (same format as the Python package).
%   tcd_save_reference(ref, '../refs/ps_D65.csv')                  % also writes ../refs/ps_D65.json
%   tcd_save_reference(ref, '../refs/ps_D65.csv', 'Spectra', true) % + ../refs/ps_D65_spectra.csv
%   The CSV holds one row per candidate: sweep parameters, label, XYZ, J'a'b', L*a*b*, sRGB, Delta E.
%   The JSON holds the label description, the stack and the system file it came from.
arguments
    ref struct
    file (1,:) char
    opts.Spectra (1,1) logical = false
end
pk = fieldnames(ref.params)';
cols = [pk, {'X', 'Y', 'Z', 'Jp', 'ap', 'bp', 'L', 'a', 'b', 'sR', 'sG', 'sB', 'dE_substrate'}];
P = zeros(numel(ref.value), numel(pk));
for i = 1:numel(pk)
    P(:, i) = ref.params.(pk{i});
end
M = [P, ref.XYZ, ref.Jab, ref.Lab, ref.sRGB, ref.dE];
folder = fileparts(file);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
fid = fopen(file, 'w');
fprintf(fid, '%s\n', strjoin(cols, ','));
fmt = [strjoin(repmat({'%.6g'}, 1, size(M, 2)), ','), '\n'];
fprintf(fid, fmt, M');
fclose(fid);

meta = struct('format', 1, 'name', ref.name, 'label', ref.label);
if isfield(ref, 'meta')
    for f = fieldnames(ref.meta)'
        meta.(f{1}) = ref.meta.(f{1});
    end
end
meta.columns = cols;
[p, n] = fileparts(file);
fid = fopen(fullfile(p, [n '.json']), 'w');
fprintf(fid, '%s\n', jsonencode(meta, 'PrettyPrint', true));
fclose(fid);

if opts.Spectra && isfield(ref, 'spectra')
    S = [(360:830)', ref.spectra'];
    fid = fopen(fullfile(p, [n '_spectra.csv']), 'w');
    fprintf(fid, 'wavelength_nm,%s\n', strjoin(arrayfun(@(i) sprintf('row%d', i - 1), 1:size(ref.spectra, 1), ...
        'UniformOutput', false), ','));
    fprintf(fid, [strjoin(repmat({'%.6g'}, 1, size(S, 2)), ','), '\n'], S');
    fclose(fid);
end
end
