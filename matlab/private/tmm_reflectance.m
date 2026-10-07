function R = tmm_reflectance(N, d, wl)
%TMM_REFLECTANCE Normal-incidence power reflectance of a planar multilayer.
%   N  : M-by-W complex refractive indices, incidence medium (air) first, substrate last
%   d  : 1-by-M thicknesses in nm (first and last entries are ignored)
%   wl : 1-by-W wavelengths in nm
%   Same formalism as TransferMatrix.m (Pettersson et al., J. Appl. Phys. 86, 487, 1999):
%   S = I_01 L_1 I_12 ... L_(M-1) I_(M-1,M),  R = |S21/S11|^2, vectorised over wavelength.
[s11, s12, s21, s22] = interface(N(1, :), N(2, :));
for j = 2:size(N, 1) - 1
    ph = 2 * pi * N(j, :) * d(j) ./ wl;
    e1 = exp(-1i * ph);
    e2 = exp(1i * ph);
    s11 = s11 .* e1;  s21 = s21 .* e1;          % S * L
    s12 = s12 .* e2;  s22 = s22 .* e2;
    [i11, i12, i21, i22] = interface(N(j, :), N(j + 1, :));
    t11 = s11 .* i11 + s12 .* i21;  t12 = s11 .* i12 + s12 .* i22;   % (S L) * I
    t21 = s21 .* i11 + s22 .* i21;  t22 = s21 .* i12 + s22 .* i22;
    s11 = t11; s12 = t12; s21 = t21; s22 = t22;
end
R = abs(s21 ./ s11) .^ 2;
end

function [a, b, c, d] = interface(n1, n2)
r = (n1 - n2) ./ (n1 + n2);
t = 2 * n1 ./ (n1 + n2);
a = 1 ./ t;  b = r ./ t;  c = r ./ t;  d = 1 ./ t;
end
