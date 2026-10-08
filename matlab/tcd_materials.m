function out = tcd_materials(material, wl)
%TCD_MATERIALS Optical constants n + ik, and the list of available materials and systems.
%
%   tcd_materials                          % print the n,k library and the bundled system files
%   names = tcd_materials()                % library names as a cell array
%   n = tcd_materials('SiO2-Franta', 360:830)   % complex index on these wavelengths (nm)
%   n = tcd_materials('../my_data/MoS2.csv', wl) % from a (wavelength, n, k) file (um or nm)
%   n = tcd_materials(1.5, wl)             % constant
%
%   The library is data/nk_library.csv (a copy of Index_of_Refraction_library.xls) plus the
%   materials added with tcd_add_material (data/materials/); case, '-', '_' and spaces are
%   ignored, so 'SiO2-Franta' and 'sio2_franta' are the same.
%   Values are interpolated and extrapolated linearly (as in the original TransferMatrix code).
%   Mixtures, Cauchy and Sellmeier models are written in system files (systems/README.md).
%
%   See also TCD_ADD_MATERIAL, TCD_LOAD_SYSTEM, TCD_TMM.
if nargin >= 1
    if nargin < 2, wl = 360:830; end
    out = load_nk(material, wl);
    return
end
root = repo_root();
lib = load_nk('list');
builtin = strcmp({lib.origin}, 'data/nk_library.csv');
if nargout > 0
    out = {lib.name};
    return
end
fprintf('n,k library, usable by name in system files:\n  built in (data/nk_library.csv): %s\n', ...
    strjoin({lib(builtin).name}, ', '));
added = lib(~builtin);
if isempty(added)
    fprintf('  added (data/materials/): none yet; add one with tcd_add_material\n\n');
else
    fprintf('  added (data/materials/):\n');
    for m = added
        src = m.source; if isempty(src), src = '(not given)'; end
        fprintf('    %-16s %.0f-%.0f nm   source: %s\n', m.name, m.wl(1), m.wl(end), src);
    end
    fprintf('\n');
end
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
