function file = tcd_add_material(name, sourceFile, opts)
%TCD_ADD_MATERIAL Add optical constants (n, k) to the library used by every system file.
%
%   tcd_add_material('MoS2', 'MoS2_refractiveindexinfo.csv', 'Source', 'Beal & Hughes 1979, refractiveindex.info')
%   tcd_add_material('PMMA', 'pmma.txt', 'Source', 'own ellipsometry, 2026', 'Notes', 'spin-coated film')
%
%   The function reads the file, checks it, and saves it as data/materials/<name>.csv. From then
%   on, every system file can use the material by name ("material": "MoS2"), in MATLAB and in
%   Python. To share the material, commit the new file and open a pull request.
%
%   Accepted files: columns wavelength, n[, k] (comma, tab, space or semicolon separated; header
%   lines are skipped), and the CSV export of refractiveindex.info (an n table followed by a k
%   table). Wavelengths in nm, or in um if every value is below 50.
%
%   Checks: n > 0, k >= 0 (n + ik, k > 0 = absorption), no repeated wavelength, at least two rows.
%   Warnings: data that do not cover 360-830 nm (the missing part is extrapolated).
%
%   Options (name-value):
%     Source (required)   where the data come from: paper, database or measurement
%     Reference ('')      DOI or URL
%     Notes ('')          anything a user should know (crystal axis, film or bulk, ...)
%     Units ('auto')      'auto' | 'nm' | 'um'
%     Replace (false)     overwrite an added material of the same name (built-in names are never replaced)
%
%   Same rules and file format as 'python python/run_tcd.py add-material'.
%   See also TCD_MATERIALS, TCD_LOAD_SYSTEM.
arguments
    name (1,:) char
    sourceFile (1,:) char
    opts.Source (1,:) char
    opts.Reference (1,:) char = ''
    opts.Notes (1,:) char = ''
    opts.Units (1,:) char {mustBeMember(opts.Units, {'auto', 'nm', 'um'})} = 'auto'
    opts.Replace (1,1) logical = false
end
if isempty(regexp(name, '^[A-Za-z][A-Za-z0-9_-]*$', 'once'))
    error('tcd:material', 'name ''%s'': use letters, digits, ''-'' and ''_'', starting with a letter (e.g. MoS2_bulk)', name);
end
if ~isfield(opts, 'Source') || isempty(strtrim(opts.Source))
    error('tcd:material', 'Give the source of the data (''Source'', ...) so others can trust it');
end
lib = load_nk('list');
hit = find(strcmp(lower(regexprep({lib.name}, '[^A-Za-z0-9]', '')), lower(regexprep(name, '[^A-Za-z0-9]', ''))), 1);
if ~isempty(hit)
    if strcmp(lib(hit).origin, 'data/nk_library.csv')
        error('tcd:material', '''%s'' is already a built-in library material (%s); choose another name', name, lib(hit).name);
    elseif ~opts.Replace
        error('tcd:material', '''%s'' already exists (%s); use ''Replace'', true to overwrite it', name, lib(hit).origin);
    end
end
[wl, n, k] = parse_nk_file(sourceFile, opts.Units);
problems = {};
if numel(wl) < 2, problems{end + 1} = 'need at least two wavelengths'; end
if any(~isfinite([wl; n; k])), problems{end + 1} = 'the table contains empty or non-numeric values'; end
if any(n <= 0), problems{end + 1} = 'n must be positive'; end
if any(k < 0), problems{end + 1} = 'k must be zero or positive (this code uses n + ik, with k > 0 for absorption)'; end
if ~isempty(problems)
    error('tcd:material', '%s: %s', sourceFile, strjoin(problems, '; '));
end
if wl(1) > 360, warning('tcd:material', 'data start at %g nm: 360-%g nm will be extrapolated', wl(1), wl(1)); end
if wl(end) < 830, warning('tcd:material', 'data end at %g nm: %g-830 nm will be extrapolated', wl(end), wl(end)); end
if wl(1) > 400 || wl(end) < 700
    warning('tcd:material', 'the data do not cover 400-700 nm: the simulated colours will be unreliable');
end
folder = fullfile(repo_root(), 'data', 'materials');
if ~isfolder(folder), mkdir(folder); end
file = fullfile(folder, [name '.csv']);
[~, srcName, srcExt] = fileparts(sourceFile);
fid = fopen(file, 'w');
fprintf(fid, '# name: %s\n# source: %s\n', name, strtrim(opts.Source));
if ~isempty(strtrim(opts.Reference)), fprintf(fid, '# reference: %s\n', strtrim(opts.Reference)); end
if ~isempty(strtrim(opts.Notes)), fprintf(fid, '# notes: %s\n', strtrim(opts.Notes)); end
fprintf(fid, '# added: %s\n# imported from: %s%s\nwavelength_nm,n,k\n', char(datetime('today', 'Format', 'yyyy-MM-dd')), srcName, srcExt);
fprintf(fid, '%.6g,%.6g,%.6g\n', [wl, n, k]');
fclose(fid);
load_nk('clear-cache');
v = load_nk(name, [450 550 650]);
fprintf('Added ''%s'' -> %s\n  450 nm: n = %.3f, k = %.3f   550 nm: n = %.3f, k = %.3f   650 nm: n = %.3f, k = %.3f\n', ...
    name, file, real(v(1)), imag(v(1)), real(v(2)), imag(v(2)), real(v(3)), imag(v(3)));
fprintf('  Use it in a system file as  "material": "%s"\n  To share it: commit data/materials/%s.csv and open a pull request.\n', name, name);
end
