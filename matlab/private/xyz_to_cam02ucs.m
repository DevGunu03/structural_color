function Jab = xyz_to_cam02ucs(XYZ)
%XYZ_TO_CAM02UCS CAM02-UCS J', a', b' of N-by-3 XYZ (0-100 scale, D65-referred).
%   CIECAM02 forward model (CIE 159:2004) followed by Luo et al. (2006) UCS, with the
%   viewing conditions of Cobeldick's sRGB_to_CAM02UCS.m. Matches colour-science.
c = colour_const();
Xw = c.XYZ_D65;  Yw = Xw(2);  L_A = c.L_A;
D = min(max(c.F * (1 - (1 / 3.6) * exp((-L_A - 42) / 92)), 0), 1);
k = 1 / (5 * L_A + 1);
F_L = 0.2 * k^4 * (5 * L_A) + 0.1 * (1 - k^4)^2 * (5 * L_A)^(1 / 3);
n = c.Y_b / Yw;
z = 1.48 + sqrt(n);
N_bb = 0.725 * (1 / n)^0.2;  N_cb = N_bb;

LMS_w = (c.M_CAT02 * Xw')';
cad = Yw * D ./ LMS_w + 1 - D;                          % chromatic adaptation factors
toHPE = c.M_HPE / c.M_CAT02;
resp = @(LMS) compress((toHPE * (LMS .* cad)')', F_L);
RGBa_w = resp(LMS_w);
A_w = (2 * RGBa_w(1) + RGBa_w(2) + RGBa_w(3) / 20 - 0.305) * N_bb;

RGBa = resp(XYZ * c.M_CAT02');
a = RGBa(:, 1) - 12 * RGBa(:, 2) / 11 + RGBa(:, 3) / 11;
b = (RGBa(:, 1) + RGBa(:, 2) - 2 * RGBa(:, 3)) / 9;
h = mod(atan2d(b, a), 360);
e_t = 0.25 * (cos(h * pi / 180 + 2) + 3.8);
A = (2 * RGBa(:, 1) + RGBa(:, 2) + RGBa(:, 3) / 20 - 0.305) * N_bb;
J = 100 * spow(A / A_w, c.c * z);
t = (50000 / 13) * c.N_c * N_cb * e_t .* sqrt(a .^ 2 + b .^ 2) ./ (RGBa(:, 1) + RGBa(:, 2) + 21 * RGBa(:, 3) / 20);
C = spow(t, 0.9) .* spow(J / 100, 0.5) * (1.64 - 0.29^n)^0.73;
M = C * F_L^0.25;

Jp = (1 + 100 * c.c1) * J ./ (1 + c.c1 * J);
Mp = log(1 + c.c2 * M) / c.c2;
Jab = [Jp, Mp .* cosd(h), Mp .* sind(h)];
end

function y = compress(x, F_L)
f = spow(F_L * x / 100, 0.42);
y = 400 * sign(x) .* f ./ (27.13 + f) + 0.1;
end

function y = spow(x, p)
y = sign(x) .* abs(x) .^ p;
end
