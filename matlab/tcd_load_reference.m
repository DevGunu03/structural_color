function ref = tcd_load_reference(file)
%TCD_LOAD_REFERENCE Read a reference CSV (+ its JSON) written by MATLAB or by the Python package.
%   ref = tcd_load_reference(fullfile('..', 'refs', 'ps_D65.csv'))
%   References written before system files existed ("system": "PS" / "MoO3") are understood too.
[p, n] = fileparts(file);
meta = jsondecode(fileread(fullfile(p, [n '.json'])));
t = readtable(file, 'VariableNamingRule', 'preserve');
names = t.Properties.VariableNames;
nPar = find(strcmp(names, 'X'), 1) - 1;
params = struct();
for i = 1:nPar
    params.(names{i}) = t{:, i};
end
if isfield(meta, 'name'), name = meta.name; else, name = meta.system; end
if isfield(meta, 'label')
    label = meta.label;
    label.classes = logical(label.classes);
elseif strcmp(name, 'MoO3')
    label = struct('name', 'thickness_nm', 'unit', 'nm', 'title', 'MoO3 thickness', 'classes', false, ...
        'class_name', '{}', 'zero_name', '', 'gap', 40);
else
    label = struct('name', 'eff_layers', 'unit', 'layers', 'title', 'Layer number', 'classes', true, ...
        'class_name', '{}L', 'zero_name', 'substrate', 'gap', 0.5);
end
ref = finish_reference(name, params, label, [t.X, t.Y, t.Z]);
ref.meta = rmfield(meta, intersect(fieldnames(meta), {'format', 'name', 'system', 'label', 'columns'}));
end
