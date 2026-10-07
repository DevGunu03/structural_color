function ref = tcd_load_reference(file)
%TCD_LOAD_REFERENCE Read a reference CSV (+ its JSON) written by MATLAB or by the Python package.
%   ref = tcd_load_reference(fullfile('..', 'refs', 'ps_D65.csv'))
[p, n] = fileparts(file);
meta = jsondecode(fileread(fullfile(p, [n '.json'])));
t = readtable(file, 'VariableNamingRule', 'preserve');
names = t.Properties.VariableNames;
nPar = find(strcmp(names, 'X'), 1) - 1;
params = struct();
for i = 1:nPar
    params.(names{i}) = t{:, i};
end
if strcmp(meta.system, 'MoO3')
    value = params.thickness_nm;
else
    value = params.eff_layers;
end
ref = finish_reference(meta.system, params, value, [t.X, t.Y, t.Z]);
ref.meta = rmfield(meta, intersect(fieldnames(meta), {'system', 'columns'}));
end
