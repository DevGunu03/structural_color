function ref = tcd_build_reference(system, opts)
%TCD_BUILD_REFERENCE Simulated colour reference for MoO3 flakes or PS-bead layers.
%
%   ref = tcd_build_reference('MoO3')                       % 0-600 nm in 1 nm steps
%   ref = tcd_build_reference('PS', 'MaxLayers', 5)          % substrate + 1-5 layers x packing
%   ref = tcd_build_reference('PS', 'Out', '../refs/ps_D65_matlab.csv')
%
%   Stacks (light from the left), normal incidence:
%     MoO3 : Air | MoO3 (t)                    | SiO2 (Oxide nm) | Si
%     PS   : Air | PS layers N..1 (top packing) | SiO2 (Oxide nm) | Si
%   PS layers are Maxwell-Garnett slabs of beads in air. Row 1 is always the bare substrate.
%   Each row's reflectance becomes XYZ (CIE 1931 2 deg, Illuminant), CAM02-UCS, CIELAB, sRGB.
%
%   Options (name-value):
%     TMax (600), Step (1)                MoO3 thickness range in nm
%     MaxLayers (5), PackingStep (0.01)   PS layer count and top-layer packing step (fraction)
%     Bead (300)                          PS bead diameter, nm
%     Model ('slab')                      'slab': every layer a Bead-thick slab at f = 0.6046
%                                         'hcp' : layers above the first at pitch Bead*sqrt(2/3), f = 0.7405
%     Oxide (100)                         SiO2 thickness, nm
%     Material                            MoO3 or PS optical constants (default 'aMoO3' / 'PS-beads')
%     Substrate ({'SiO2-Franta','Si-Franta'})
%     Illuminant ('D65')                  'D65' | 'A' | 'D50' | CSV of a measured lamp spectrum
%     Out ('')                            if given, also save CSV + JSON (readable by the Python code)
%
%   See also TCD_MAP_IMAGE, TCD_LOAD_REFERENCE.
arguments
    system (1,:) char {mustBeMember(system, {'MoO3', 'PS'})}
    opts.TMax (1,1) double = 600
    opts.Step (1,1) double = 1
    opts.MaxLayers (1,1) double = 5
    opts.PackingStep (1,1) double = 0.01
    opts.Bead (1,1) double = 300
    opts.Model (1,:) char {mustBeMember(opts.Model, {'slab', 'hcp'})} = 'slab'
    opts.Oxide (1,1) double = 100
    opts.Material (1,:) char = ''
    opts.Substrate cell = {'SiO2-Franta', 'Si-Franta'}
    opts.Illuminant (1,:) char = 'D65'
    opts.Out (1,:) char = ''
end
c = colour_const();
wl = c.wl;
nAir = load_nk('Air', wl);
nOx = load_nk(opts.Substrate{1}, wl);
nSi = load_nk(opts.Substrate{2}, wl);

if strcmp(system, 'MoO3')
    mat = opts.Material; if isempty(mat), mat = 'aMoO3'; end
    nF = load_nk(mat, wl);
    t = (0:opts.Step:opts.TMax)';
    R = zeros(numel(t), numel(wl));
    for i = 1:numel(t)
        R(i, :) = tmm_reflectance([nAir; nF; nOx; nSi], [0 t(i) opts.Oxide 0], wl);
    end
    params = struct('thickness_nm', t);
    value = t;
    stack = sprintf('Air | %s (t) | %s %g nm | %s', mat, opts.Substrate{1}, opts.Oxide, opts.Substrate{2});
else
    mat = opts.Material; if isempty(mat), mat = 'PS-beads'; end
    nPS = load_nk(mat, wl);
    packs = round((opts.PackingStep:opts.PackingStep:1 + opts.PackingStep / 2)', 6);
    packs = packs(packs <= 1 + 1e-9);
    [P, L] = ndgrid(packs, 1:opts.MaxLayers);
    layers = [0; L(:)];
    packing = [0; P(:)];
    R = zeros(numel(layers), numel(wl));
    for i = 1:numel(layers)
        [f, d] = ps_stack(layers(i), packing(i), opts.Bead, opts.Model);
        N = nAir;
        for j = 1:numel(f)
            N = [N; maxwell_garnett(nAir, nPS, f(j))]; %#ok<AGROW>
        end
        R(i, :) = tmm_reflectance([N; nOx; nSi], [0 d opts.Oxide 0], wl);
    end
    eff = layers - 1 + packing;
    eff(layers == 0) = 0;
    params = struct('layers', layers, 'packing', packing, 'eff_layers', eff);
    value = eff;
    stack = sprintf('Air | %s layers (d=%g nm, model=%s) | %s %g nm | %s', mat, opts.Bead, opts.Model, ...
        opts.Substrate{1}, opts.Oxide, opts.Substrate{2});
end

ref = finish_reference(system, params, value, refl_to_xyz(R, opts.Illuminant));
ref.spectra = R;
ref.meta = struct('stack', stack, 'illuminant', opts.Illuminant, 'NA', 0, 'created_by', 'MATLAB tcd_build_reference');
if ~isempty(opts.Out)
    tcd_save_reference(ref, opts.Out);
end
end

function [f, d] = ps_stack(n, p, bead, model)
% Fill fractions and thicknesses of the PS slabs, top layer first.
f = zeros(1, n); d = zeros(1, n);
for k = n:-1:1                     % k = layer index counted from the substrate
    if strcmp(model, 'slab') || k == 1
        fk = pi / (3 * sqrt(3));  dk = bead;
    else
        fk = pi / (3 * sqrt(2));  dk = bead * sqrt(2 / 3);
    end
    if k == n, fk = fk * p; end
    f(n - k + 1) = fk;  d(n - k + 1) = dk;
end
end
