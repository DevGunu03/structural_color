function ok = run_tests()
%RUN_TESTS Self-contained checks of the MATLAB TCD code (no lab images needed).
%   From the matlab/ folder:  addpath tests; run_tests
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fileparts(here));
fx = fullfile(root, 'tests', 'fixtures');
wl = 360:830;
fails = 0;

% 1. transfer matrix against the original TransferMatrix_multiple/_packing.m output
nk = @(name) tcd_materials(name, wl);
d = readmatrix(fullfile(fx, 'matlab_MoO3_100nm.txt'));
R = tcd_tmm([nk('Air'); nk('aMoO3'); nk('SiO2-Franta'); nk('Si-Franta')], [0 100 100 0], wl);
fails = fails + check('TMM vs TransferMatrix_multiple.m, MoO3 100 nm', max(abs(R - d(:, 2)')) < 1e-10);
d = readmatrix(fullfile(fx, 'matlab_PS_solid_film_600nm.txt'));
R = tcd_tmm([nk('Air'); nk('PS-beads'); nk('SiO2-Franta'); nk('Si-Franta')], [0 600 100 0], wl);
fails = fails + check('TMM vs TransferMatrix_multiple.m, solid PS 600 nm', max(abs(R - d(:, 2)')) < 1e-10);
S = tcd_spectrum('ps', struct('layers', 2, 'packing', 0.5));
d = readmatrix(fullfile(fx, 'matlab_PS_bilayer_top50pct.txt'));
fails = fails + check('systems/ps.jsonc vs TransferMatrix_packing.m, PS bilayer 50 %', max(abs(S.R - d(:, 2)')) < 1e-5);
N = [nk('Air'); nk('aMoO3'); nk('SiO2-Franta'); nk('Si-Franta')];
d = readmatrix(fullfile(fx, 'python_MoO3_100nm_NA0.5.txt'));
fails = fails + check('NA 0.5 cone average equals the Python code', ...
    max(abs(tcd_tmm(N, [0 100 100 0], wl, 'NA', 0.5) - d(:, 2)')) < 1e-10);
fails = fails + check('NA -> 0 limit of the cone average', ...
    max(abs(tcd_tmm(N, [0 100 100 0], wl, 'NA', 0.02) - tcd_tmm(N, [0 100 100 0], wl))) < 1e-3);
fails = fails + check('zero-thickness layer changes nothing', ...
    max(abs(tcd_tmm(N, [0 0 100 0], wl) - tcd_tmm(N([1 3 4], :), [0 100 0], wl))) < 1e-12);

% 2. colour conversions
jab = tcd_colour('srgb2cam02ucs', [0.34712092 0.44539356 0.59242291]);
fails = fails + check('CAM02-UCS vs CIECAM02 toolbox', max(abs(jab - [49.19585505 -5.70443096 -17.40028325])) < 0.01);
moo3 = tcd_load_reference(fullfile(root, 'refs', 'moo3_D65.csv'));
fails = fails + check('CIELAB L* of bare substrate', moo3.Lab(1, 1) > 40 && moo3.Lab(1, 1) < 60);

% 3. expression language (same rules as python/tcd/expr.py), through a system's constants
exprs = {'1 + 2 * 3', 7; '-2^2', -4; '2^-1', 0.5; '2**3', 8; '(1 + 2) * 3', 9; 'max(1, min(5, 3))', 3; ...
    'round(2.5) + round(-2.5)', 1; '3 >= 2', 1; '2 == 3', 0; 'sqrt(16) / 2e0', 2; 'pi', pi; ...
    '1 - 2 - 3', -4; '8 / 4 / 2', 1; '2^3^2', 512; '.5 + 1e-3', 0.501};
def = struct('name', 'expr_test', 'constants', struct(), 'sweep', struct('t', 0), 'layers', {{}}, 'substrate', 'Si-Franta');
for i = 1:size(exprs, 1), def.constants.(sprintf('c%d', i)) = exprs{i, 1}; end
s = tcd_load_system(def);
got = cellfun(@(f) s.constants.(f), fieldnames(s.constants));
fails = fails + check('expressions: precedence, functions, comparisons', max(abs(got - [exprs{:, 2}]')) < 1e-12);
bad = {'t + 1', 'unknown name'; '2 +', 'end of expression'; 'max(1)', 'argument'; 'system(''dir'')', ''};
okExpr = true;
for i = 1:size(bad, 1)
    d2 = def; d2.constants = struct('c', bad{i, 1});
    okExpr = okExpr && throws(@() tcd_load_system(d2), bad{i, 2});
end
fails = fails + check('expressions: unknown names and code are refused', okExpr);

% 4. system files
files = dir(fullfile(root, 'systems', '*.jsonc'));
nLoaded = 0;
for k = 1:numel(files)
    try
        tcd_load_system(fullfile(files(k).folder, files(k).name));
        nLoaded = nLoaded + 1;
    catch err
        fprintf('    %s: %s\n', files(k).name, err.message);
    end
end
fails = fails + check(sprintf('all %d bundled system files load', numel(files)), nLoaded == numel(files));
for f = {'moo3', 'ps', 'ps_hcp', 'sio2_on_si', 'graphene', 'bragg_mirror'}
    py = tcd_load_reference(fullfile(root, 'refs', [f{1} '_D65.csv']));
    mm = tcd_build_reference(f{1});
    same = all(cellfun(@(c) max(abs(mm.params.(c) - py.params.(c))) < 1e-6, fieldnames(py.params)));
    fails = fails + check(sprintf('MATLAB reference equals Python refs/%s_D65.csv', f{1}), ...
        same && max(abs(mm.Jab - py.Jab), [], 'all') < 1e-3 && strcmp(mm.meta.stack, py.meta.stack));
end
s = tcd_load_system('ps', 'Set', {'oxide_nm', 285; 'max_layers', 2});
fails = fails + check('Set overrides constants', s.constants.oxide_nm == 285 && max(s.rows.layers) == 2);
S = tcd_spectrum('bragg_mirror', struct('lambda0', 550), 'Set', {'pairs', 2});
hi = nk('TiO2-Franta');  lo = nk('SiO2-Franta');
R = tcd_tmm([nk('Air'); hi; lo; hi; lo; nk('Si-Franta')], [0 550/4/2.35 550/4/1.46 550/4/2.35 550/4/1.46 0], wl);
fails = fails + check('a repeated group equals the layers written out', max(abs(S.R - R)) < 1e-12);
base = tcd_load_system('moo3');
base = base.definition;
typo = base;        typo.thicknes = 1;
unknownMat = base;  unknownMat.substrate = 'Unobtainium';
overfull = base;
overfull.materials = struct('m', struct('ema', 'bruggeman', 'host', 'Air', 'inclusion', 'aMoO3', ...
    'fraction', 'thickness_nm / 100'));
overfull.layers = {struct('material', 'm', 'thickness', 5)};
badExpr = base;     badExpr.layers = {struct('material', 'aMoO3', 'thickness', 'thickness_nm +')};
mistakes = {typo, 'unknown key'; unknownMat, 'n,k library'; overfull, 'within [0, 1]'; badExpr, 'thickness'};
okSys = all(cellfun(@(d, t) throws(@() tcd_load_system(d), t), mistakes(:, 1), mistakes(:, 2)));
okSys = okSys && throws(@() tcd_load_system('moo3', 'Set', {'oxide', 285}), 'not a constant');
fails = fails + check('system file mistakes give clear errors', okSys);

% 4b. reading n,k files and adding materials to the library
tmp = tempname;  mkdir(tmp);
write_lines(fullfile(tmp, 'ri.csv'), {'wl,n', '0.30,2.0', '0.50,2.1', '0.90,2.2', 'wl,k', '0.30,0.3', '0.90,0.1'});
write_lines(fullfile(tmp, 'negk.csv'), {'400,2,-0.1', '800,2,0'});
nri = tcd_materials(fullfile(tmp, 'ri.csv'), [300 500 900]);
fails = fails + check('n,k files: refractiveindex.info n + k blocks in um', ...
    max(abs(nri - ([2 2.1 2.2] + 1i * [0.3 0.7 / 3 0.1]))) < 1e-9);
name = 'ZZ_test_material_m';
target = fullfile(root, 'data', 'materials', [name '.csv']);
try
    tcd_add_material(name, fullfile(tmp, 'ri.csv'), 'Source', 'synthetic test data');
    okMat = isfile(target) && any(strcmp(tcd_materials(), name));
    def = tcd_load_system('sio2_on_si');  def = def.definition;
    def.layers = {struct('name', 'film', 'material', name, 'thickness', 'oxide_nm')};
    s = tcd_load_system(def);
    okMat = okMat && s.nRows == 1001;
    okMat = okMat && throws(@() tcd_add_material(name, fullfile(tmp, 'ri.csv'), 'Source', 'x'), 'already exists');
    okMat = okMat && throws(@() tcd_add_material('SiO2-Franta', fullfile(tmp, 'ri.csv'), 'Source', 'x'), 'built-in');
    okMat = okMat && throws(@() tcd_add_material('bad name', fullfile(tmp, 'ri.csv'), 'Source', 'x'), 'letters');
    okMat = okMat && throws(@() tcd_add_material('ZZ_other', fullfile(tmp, 'negk.csv'), 'Source', 'x'), 'k must be');
catch err
    fprintf('    %s\n', err.message);
    okMat = false;
end
if isfile(target), delete(target); end
tcd_materials('Air', 500);  clear load_nk;     % forget the deleted material
fails = fails + check('tcd_add_material: saved to data/materials, usable by name, mistakes refused', okMat);
rmdir(tmp, 's');

% 5. physics sanity: monolayer graphene contrast peaks at the oxide's reflectance minimum
S = tcd_spectrum('graphene', struct('layers', [0 1]), 'Set', {'oxide_nm', 300});
vis = wl >= 450 & wl <= 700;
C = (S.R(1, :) - S.R(2, :)) ./ S.R(1, :);
[cmax, i1] = max(C(vis));  [~, i2] = min(S.R(1, vis));  wv = wl(vis);
fails = fails + check(sprintf('graphene on 300 nm SiO2: contrast %.3f at %d nm (R min at %d nm)', cmax, wv(i1), wv(i2)), ...
    cmax > 0.08 && cmax < 0.16 && abs(wv(i1) - wv(i2)) < 15);

% 6. reference files: round trip, and references written before system files
tmp = tempname;  mkdir(tmp);
g = tcd_build_reference('graphene');
tcd_save_reference(g, fullfile(tmp, 'g.csv'));
back = tcd_load_reference(fullfile(tmp, 'g.csv'));
fails = fails + check('reference save/load keeps the label description', isequal(back.label, g.label));
meta = jsondecode(fileread(fullfile(root, 'refs', 'ps_D65.json')));
legacy = struct('stack', meta.stack, 'illuminant', 'D65', 'NA', 0, 'system', 'PS');
fid = fopen(fullfile(tmp, 'old.json'), 'w'); fprintf(fid, '%s', jsonencode(legacy)); fclose(fid);
copyfile(fullfile(root, 'refs', 'ps_D65.csv'), fullfile(tmp, 'old.csv'));
old = tcd_load_reference(fullfile(tmp, 'old.csv'));
fails = fails + check('old-style PS reference still loads', old.label.classes && strcmp(old.label.name, 'eff_layers'));
rmdir(tmp, 's');

% 7. synthetic micrographs: simulated colours, unknown camera white balance, noise
XYZ2RGB = inv([0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505]);
gains = [0.55 0.80 1.25];
ps = tcd_load_reference(fullfile(root, 'refs', 'ps_D65.csv'));
rows = [1, arrayfun(@(k) find(ps.params.layers == k & abs(ps.params.packing - 1) < 1e-9), 1:3)];
file = synth(ps.XYZ(rows, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, ps, 'Substrate', [11 11 30 30], 'ShowFigure', false);
delete(file);
got = arrayfun(@(b) mode(r.classes(21:60, (b - 1) * 60 + 16:b * 60 - 15), 'all'), 1:4);
fails = fails + check(sprintf('synthetic PS image: layers %s', mat2str(got)), isequal(got, [0 1 2 3]));

t = [0 150 280 420];
file = synth(moo3.XYZ(t + 1, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, moo3, 'Substrate', [11 11 30 30], 'ShowFigure', false);
delete(file);
got = arrayfun(@(b) median(r.value(21:60, (b - 1) * 60 + 16:b * 60 - 15), 'all'), 1:4);
fails = fails + check(sprintf('synthetic MoO3 image: %s nm', mat2str(got)), all(abs(got - t) <= 3));

gr = tcd_load_reference(fullfile(root, 'refs', 'graphene_D65.csv'));
truth = [0 1 3 6];
file = synth(gr.XYZ(truth + 1, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, gr, 'Substrate', [11 11 30 30], 'ShowFigure', false);
delete(file);
got = arrayfun(@(b) mode(r.classes(21:60, (b - 1) * 60 + 16:b * 60 - 15), 'all'), 1:4);
fails = fails + check(sprintf('synthetic graphene image: layers %s', mat2str(got)), isequal(got, truth));

% checks and reliability on a synthetic MoO3 image (the model is exact here)
file = synth(moo3.XYZ(t + 1, :) / 100 * XYZ2RGB', gains, tempname);
r = tcd_map_image(file, moo3, 'Substrate', [11 11 30 30], 'ShowFigure', false, ...
    'Checks', {'150nm', [76 26 20 30]; '400nm', [136 26 20 30]});
delete(file);
c = r.checks;
fails = fails + check(sprintf('checks: right one agrees (%.0f %%), wrong one is flagged (%.0f %%), model error from checks %.2f dE', ...
    100 * c.within_tolerance(1), 100 * c.within_tolerance(2), r.sigma.model), ...
    c.within_tolerance(1) > 0.9 && c.within_tolerance(2) < 0.1 && strcmp(r.sigma.model_source, 'checks') && r.sigma.model < 1.5);
band = r.reliability(21:60, 76:105);
rel = r.reliability(~isnan(r.reliability));
fails = fails + check(sprintf('reliability in [0, 1], high on the exact synthetic 150 nm band (%.2f)', median(band(:), 'omitnan')), ...
    min(rel) >= 0 && max(rel) <= 1 && median(band(:), 'omitnan') > 0.5);

% 8. camera gamma from a synthetic exposure series
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

function tf = throws(fn, text)
% True if fn() raises an error whose message contains TEXT.
try
    fn();
    tf = false;
catch err
    tf = contains(lower(err.message), lower(text));
    if ~tf, fprintf('    unexpected message: %s\n', err.message); end
end
end

function write_lines(file, lines)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
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
