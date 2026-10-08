function out = tcd_materials(material, wl)
%TCD_MATERIALS Optical constants n + ik, and the list of available materials and systems.
%
%   tcd_materials                          % print the n,k library and the bundled system files
%   names = tcd_materials()                % library names as a cell array
%   n = tcd_materials('SiO2-Franta', 360:830)   % complex index on these wavelengths (nm)
%   n = tcd_materials('../my_data/MoS2.csv', wl) % from a (wavelength, n, k) file (um or nm)
%   n = tcd_materials(1.5, wl)             % constant
%
%   Library names come from data/nk_library.csv (a copy of Index_of_Refraction_library.xls);
%   case, '-', '_' and spaces are ignored, so 'SiO2-Franta' and 'sio2_franta' are the same.
%   Values are interpolated and extrapolated linearly (as in the original TransferMatrix code).
%   Mixtures, Cauchy and Sellmeier models are written in system files (systems/README.md).
%
%   See also TCD_LOAD_SYSTEM, TCD_TMM.
if nargin >= 1
    if nargin < 2, wl = 360:830; end
    out = load_nk(material, wl);
    return
end
root = repo_root();
fid = fopen(fullfile(root, 'data', 'nk_library.csv'));
header = fgetl(fid);
fclose(fid);
cols = strsplit(header, ',');
names = regexprep(cols(endsWith(cols, '_n')), '_n$', '');
[~, o] = sort(lower(names));
names = names(o);
if nargout > 0
    out = names;
    return
end
fprintf('n,k library (data/nk_library.csv), usable by name in system files:\n  %s\n\n', strjoin(names, ', '));
fprintf('Bundled systems (systems/), usable as tcd_build_reference(''<name>''):\n');
d = dir(fullfile(root, 'systems', '*.jsonc'));
for k = 1:numel(d)
    try
        def = read_jsonc(fullfile(d(k).folder, d(k).name));
        desc = '';
        if isfield(def, 'description'), desc = def.description; end
    catch err
        desc = ['(cannot read: ' err.message ')'];
    end
    fprintf('  %-16s %s\n', regexprep(d(k).name, '\.jsonc$', ''), desc);
end
end
