function s = system_describe(sys)
%SYSTEM_DESCRIBE One-line stack, e.g. 'Air | flake: aMoO3, thickness_nm nm | oxide: SiO2-Franta, 100 nm | Si-Franta'.
%   Same text as System.describe in python/tcd/system.py.
parts = [{material_text(sys, sys.ambient)}, walk(sys, sys.layers), {material_text(sys, sys.substrate)}];
s = strjoin(parts, ' | ');
end

function parts = walk(sys, layers)
parts = {};
for j = 1:numel(layers)
    item = layers{j};
    rep = 1;
    if isfield(item, 'repeat'), rep = item.repeat; end
    if (isnumeric(rep) && rep == 1) || (ischar(rep) && strcmp(strtrim(rep), '1'))
        repTxt = '';
    else
        repTxt = sprintf(' (repeat %s)', value_text(sys, rep));
    end
    if isfield(item, 'layers')
        parts{end + 1} = sprintf('[%s]%s', strjoin(walk(sys, item.layers), ' | '), repTxt); %#ok<AGROW>
    else
        name = '';
        if isfield(item, 'name'), name = [char(item.name) ': ']; end
        t = value_text(sys, item.thickness);
        if any(t == ' '), t = ['(' t ') nm']; else, t = [t ' nm']; end
        parts{end + 1} = sprintf('%s%s, %s%s', name, material_text(sys, item.material), t, repTxt); %#ok<AGROW>
    end
end
end

function t = value_text(sys, v)
if ischar(v) || isstring(v)
    try
        v = expr_eval(char(v), sys.constants);   % only constants: show the number
    catch
        t = char(v);
        return
    end
end
t = sprintf('%.4g', v);
end

function t = material_text(sys, spec)
if isnumeric(spec)
    t = sprintf('n=%g', spec);
elseif ischar(spec) || isstring(spec)
    t = char(spec);
elseif isfield(spec, 'library')
    t = char(spec.library);
elseif isfield(spec, 'file')
    [~, n, e] = fileparts(char(spec.file));
    t = [n e];
elseif isfield(spec, 'n')
    t = ['n=' num2str(spec.n)];
    if isfield(spec, 'k') && ~(isnumeric(spec.k) && spec.k == 0), t = [t '+' num2str(spec.k) 'i']; end
elseif isfield(spec, 'cauchy')
    t = 'Cauchy';
elseif isfield(spec, 'sellmeier')
    t = 'Sellmeier';
else
    t = sprintf('%s(%s in %s, f=%s)', spec.ema, material_text(sys, spec.inclusion), ...
        material_text(sys, spec.host), value_text(sys, spec.fraction));
end
end
