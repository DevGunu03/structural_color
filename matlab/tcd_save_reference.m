function tcd_save_reference(ref, file)
%TCD_SAVE_REFERENCE Write a reference as CSV + JSON (same format as the Python package).
%   tcd_save_reference(ref, 'refs/ps_D65.csv') also writes refs/ps_D65.json.
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
meta = ref.meta;
meta.system = ref.system;
meta.columns = cols;
[p, n] = fileparts(file);
fid = fopen(fullfile(p, [n '.json']), 'w');
fprintf(fid, '%s', jsonencode(meta, 'PrettyPrint', true));
fclose(fid);
end
