function ref = tcd_build_reference(system, opts)
%TCD_BUILD_REFERENCE Simulated colour reference of a system file: every candidate structure's colour.
%
%   ref = tcd_build_reference('moo3')                                  % ../systems/moo3.jsonc
%   ref = tcd_build_reference('ps', 'Set', {'oxide_nm', 285}, 'Out', 'auto')
%   ref = tcd_build_reference('../systems/my_sample.jsonc', 'Out', '../refs/my_sample_D65.csv')
%   ref = tcd_build_reference(sys)                                     % from tcd_load_system
%
%   For every candidate structure of the system (each row of its sweep) this builds the layer
%   stack, computes its reflectance spectrum with tcd_tmm, and turns it into a colour:
%   XYZ (CIE 1931 2 deg observer under the illuminant, then CAT02 to D65 if needed),
%   CAM02-UCS J'a'b', CIELAB and sRGB. Row 1 is the calibration structure (normally the bare
%   substrate); the 'label' column is what tcd_map_image reports for each pixel.
%
%   Options (name-value):
%     Set ({})          change constants of the system file: {'oxide_nm', 285; 'bead_nm', 500}
%     Illuminant ('')   override optics.illuminant: 'D65' | 'A' | 'D50' | CSV of a lamp spectrum
%     NA ([])           override optics.na (objective numerical aperture; 0 = normal incidence)
%     Out ('')          save CSV + JSON (readable by the Python code) + the reference sheet PNG;
%                       'auto' = ../refs/<name>[_<changes>]_<illuminant>.csv
%     Spectra (false)   also save the reflectance spectra (<Out>_spectra.csv)
%
%   ref fields: name, params (sweep columns + label), label, value, XYZ, Jab, Lab, sRGB, dE
%   (Delta E from row 1), spectra (rows x 471, 360:830 nm), meta.
%
%   The old calls tcd_build_reference('MoO3') and ('PS') still work: they load the bundled
%   systems/moo3.jsonc and systems/ps.jsonc.
%
%   See also TCD_LOAD_SYSTEM, TCD_SPECTRUM, TCD_MAP_IMAGE, TCD_PLOT_REFERENCE, TCD_TMM.
arguments
    system
    opts.Set = {}
    opts.Illuminant (1,:) char = ''
    opts.NA double = []
    opts.Out (1,:) char = ''
    opts.Spectra (1,1) logical = false
end
if isstruct(system) && isfield(system, 'rows')
    sys = system;
else
    sys = tcd_load_system(system, 'Set', opts.Set);
end
illuminant = opts.Illuminant;
if isempty(illuminant), illuminant = sys.optics.illuminant; end
na = opts.NA;
if isempty(na), na = sys.optics.na; end

wl = 360:830;
R = zeros(sys.nRows, numel(wl));
for i = 1:sys.nRows
    [N, d] = system_stack(sys, i, wl);
    R(i, :) = tcd_tmm(N, d, wl, 'NA', na, 'Angles', sys.optics.angles, 'Weighting', sys.optics.na_weighting);
end

ref = finish_reference(sys.name, sys.rows, sys.label, refl_to_xyz(R, illuminant));
ref.spectra = R;
sysFile = sys.file;
root = repo_root();
if startsWith(sysFile, root)
    sysFile = strrep(sysFile(numel(root) + 2:end), '\', '/');
end
ref.meta = struct('description', sys.description, 'stack', system_describe(sys), 'illuminant', illuminant, ...
    'NA', na, 'na_weighting', sys.optics.na_weighting, 'constants', sys.constants, 'system_file', sysFile, ...
    'created_by', 'MATLAB tcd_build_reference', 'definition', sys.definition);

if ~isempty(opts.Out)
    out = opts.Out;
    if strcmpi(out, 'auto'), out = default_out(sys, opts.Set, illuminant, na); end
    tcd_save_reference(ref, out, 'Spectra', opts.Spectra);
    [p, n] = fileparts(out);
    tcd_plot_reference(ref, fullfile(p, [n '.png']), false);
    fprintf('%s: %d candidates -> %s (+ .json, reference sheet .png)\n', sys.name, sys.nRows, out);
end
end

function out = default_out(sys, set, illuminant, na)
tag = '';
if iscell(set) && ~isempty(set)
    if size(set, 2) == 2 && ~isvector(set), pairs = set; else, pairs = reshape(set, 2, [])'; end
    for i = 1:size(pairs, 1)
        v = pairs{i, 2};
        if isnumeric(v), v = sprintf('%g', v); end
        tag = [tag '_' char(pairs{i, 1}) char(v)]; %#ok<AGROW>
    end
elseif isstruct(set)
    for f = fieldnames(set)'
        v = set.(f{1});
        if isnumeric(v), v = sprintf('%g', v); end
        tag = [tag '_' f{1} char(v)]; %#ok<AGROW>
    end
end
if na > 0, tag = sprintf('%s_NA%g', tag, na); end
[~, ill] = fileparts(illuminant);
out = fullfile(repo_root(), 'refs', regexprep([sys.name tag '_' ill '.csv'], '[^\w.-]', '_'));
end
