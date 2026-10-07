function c = colour_const()
%COLOUR_CONST Constants shared by the colour conversions (identical to the Python package).
persistent cc
if isempty(cc)
    cc.wl = 360:830;
    cmf = readmatrix(fullfile(repo_root(), 'data', 'cie1931_2deg.csv'));
    cc.cmf = cmf(:, 2:4);                       % 471-by-3, CIE 1931 2 degree
    ill = readtable(fullfile(repo_root(), 'data', 'cie_illuminants.csv'), 'VariableNamingRule', 'preserve');
    cc.illum = ill;
    cc.XYZ_D65 = [95.045592705167, 100, 108.905775075988];   % from xy = (0.3127, 0.3290)
    cc.RGB2XYZ = [0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505];  % IEC 61966-2-1
    cc.XYZ2RGB = inv(cc.RGB2XYZ);
    cc.M_CAT02 = [0.7328 0.4296 -0.1624; -0.7036 1.6975 0.0061; 0.0030 0.0136 0.9834];
    cc.M_HPE = [0.38971 0.68898 -0.07868; -0.22981 1.18340 0.04641; 0 0 1];
    % CAM02-UCS viewing conditions of sRGB_to_CAM02UCS.m (Cobeldick): L_A = 64/(5 pi), Y_b = 20, average
    cc.L_A = 64 / pi / 5;  cc.Y_b = 20;  cc.F = 1.0;  cc.c = 0.69;  cc.N_c = 1.0;
    cc.c1 = 0.007;  cc.c2 = 0.0228;
end
c = cc;
end
