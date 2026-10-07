function ok = run_tests()
%RUN_TESTS Self-contained checks of the MATLAB TCD code (no lab images needed).
%   From the matlab/ folder:  addpath tests; run_tests
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fileparts(here));
fx = fullfile(root, 'tests', 'fixtures');
fails = 0;

% 1. transfer matrix against the original TransferMatrix_multiple/_packing.m output
m = tcd_build_reference('MoO3', 'TMax', 100, 'Step', 100);
d = readmatrix(fullfile(fx, 'matlab_MoO3_100nm.txt'));
fails = fails + check('TMM vs TransferMatrix.m, MoO3 100 nm', max(abs(m.spectra(2, :) - d(:, 2)')) < 1e-10);
p = tcd_build_reference('PS', 'MaxLayers', 2, 'PackingStep', 0.5);
row = find(p.params.layers == 2 & abs(p.params.packing - 0.5) < 1e-9);
d = readmatrix(fullfile(fx, 'matlab_PS_bilayer_top50pct.txt'));
fails = fails + check('TMM vs TransferMatrix_packing.m, PS bilayer 50 %', max(abs(p.spectra(row, :) - d(:, 2)')) < 1e-5);

% 2. CAM02-UCS against Cobeldick's sRGB_to_CAM02UCS.m (value computed with the toolbox)
jab = tcd_colour('srgb2cam02ucs', [0.34712092 0.44539356 0.59242291]);
fails = fails + check('CAM02-UCS vs CIECAM02 toolbox', max(abs(jab - [49.19585505 -5.70443096 -17.40028325])) < 0.01);

% 3. CIELAB of bare 100 nm SiO2/Si is physical (the old notebooks gave L* ~ 1)
fails = fails + check('CIELAB L* of bare substrate', m.Lab(1, 1) > 40 && m.Lab(1, 1) < 60);

% 4. MATLAB references equal the Python-built ones in refs/
for f = {'moo3_D65', 'ps_D65'}
    py = tcd_load_reference(fullfile(root, 'refs', [f{1} '.csv']));
    sys = 'PS'; if startsWith(f{1}, 'moo3'), sys = 'MoO3'; end
    mm = tcd_build_reference(sys);
    fails = fails + check(['reference matches Python: ' f{1}], max(abs(mm.Jab - py.Jab), [], 'all') < 1e-3);
end

% 5-6. synthetic micrographs: simulated colours, unknown camera white balance, noise
XYZ2RGB = inv([0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505]);
gains = [0.55 0.80 1.25];
ps = tcd_load_reference(fullfile(root, 'refs', 'ps_D65.csv'));
truth = [0 1 2 3];
rows = [1, arrayfun(@(k) find(ps.params.layers == k & abs(ps.params.packing - 1) < 1e-9), 1:3)];
file = synth(ps.XYZ(rows, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, ps, 'Substrate', [11 11 30 30], 'ShowFigure', false);
delete(file);
got = arrayfun(@(b) mode(r.layers(21:60, (b - 1) * 60 + 16:b * 60 - 15), 'all'), 1:4);
fails = fails + check(sprintf('synthetic PS image: layers %s', mat2str(got)), isequal(got, truth));

mo = tcd_load_reference(fullfile(root, 'refs', 'moo3_D65.csv'));
t = [0 150 280 420];
file = synth(mo.XYZ(t + 1, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, mo, 'Substrate', [11 11 30 30], 'ShowFigure', false);
delete(file);
got = arrayfun(@(b) median(r.value(21:60, (b - 1) * 60 + 16:b * 60 - 15), 'all'), 1:4);
fails = fails + check(sprintf('synthetic MoO3 image: %s nm', mat2str(got)), all(abs(got - t) <= 3));

% 7. camera gamma from a synthetic exposure series
tt = [0.25 0.5 1 1.5 2];
for enc = {'linear', 'srgb'}
    files = cell(1, numel(tt));
    for i = 1:numel(tt)
        v = 0.3 * tt(i) * ones(60, 80, 3);
        if strcmp(enc{1}, 'srgb'), v = (v > 0.0031308) .* (1.055 * v .^ (1 / 2.4) - 0.055) + (v <= 0.0031308) .* 12.92 .* v; end
        files{i} = [tempname '.png'];
        imwrite(uint8(round(255 * min(v, 1))), files{i});
    end
    g = mean(tcd_measure_gamma(files, tt));
    delete(files{:});
    expect = 1; if strcmp(enc{1}, 'srgb'), expect = 2.2; end
    fails = fails + check(sprintf('gamma of a %s camera = %.2f', enc{1}, g), abs(g - expect) < 0.3);
end

ok = fails == 0;
if ok, fprintf('\nall checks passed\n'); else, fprintf('\n%d check(s) failed\n', fails); end
end

function f = check(name, pass)
if pass, fprintf('[PASS] %s\n', name); else, fprintf('[FAIL] %s\n', name); end
f = ~pass;
end

function file = synth(lin, gains, stem)
% Four 60-px bands of the given linear colours, as a linear camera with these channel gains sees them.
img = zeros(80, 240, 3);
rng(0);
for b = 1:4
    img(:, (b - 1) * 60 + 1:b * 60, :) = repmat(reshape(lin(b, :) .* gains, 1, 1, 3), 80, 60);
end
img = min(max(img + 0.004 * randn(size(img)), 0), 1);
file = [stem '.png'];
imwrite(uint8(round(255 * img)), file);
end
