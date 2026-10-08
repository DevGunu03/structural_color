function R = tcd_tmm(N, d, wl, opts)
%TCD_TMM Reflectance of a planar multilayer by the transfer-matrix method.
%
%   R = tcd_tmm(N, d, wl)                    normal incidence (what a low-NA objective sees)
%   R = tcd_tmm(N, d, wl, 'NA', 0.5)         unpolarised, averaged over the objective's light cone
%   R = tcd_tmm(N, d, wl, 'Theta', 0.3, 'Polarisation', 'p')   one angle (rad) and polarisation
%
%   N   M-by-W complex refractive indices n + ik (k > 0 = absorbing), one row per medium, ordered
%       from the side the light comes from: row 1 = ambient (air), row M = substrate. Both
%       outer media are semi-infinite.
%   d   1-by-M thicknesses in nm. d(1) and d(M) are ignored (semi-infinite media).
%   wl  1-by-W vacuum wavelengths in nm.
%   R   1-by-W power reflectance, between 0 and 1.
%
%   Example (100 nm MoO3 on 100 nm SiO2 on Si):
%       wl = 360:830;
%       N = [tcd_materials('Air', wl); tcd_materials('aMoO3', wl);
%            tcd_materials('SiO2-Franta', wl); tcd_materials('Si-Franta', wl)];
%       R = tcd_tmm(N, [0 100 100 0], wl);   plot(wl, R)
%   You rarely call this directly: tcd_build_reference and tcd_spectrum build N and d from a
%   system file (systems/*.jsonc) and call it for every candidate structure.
%
%   Options (name-value):
%     NA (0)                 numerical aperture of the objective (in air). 0 = normal incidence.
%     Angles (24)            Gauss-Legendre points over the cone [0, asin(NA)]
%     Weighting ('uniform')  'uniform' (filled back aperture) or 'gaussian' (w = exp(-(theta/theta_max)^2))
%     Theta ([])             a single angle of incidence in the ambient (rad); overrides NA
%     Polarisation ('s')     's' or 'p', used with Theta
%
%   ----------------------------------------------------------------------------------------
%   The method (Pettersson, Roman & Inganas, J. Appl. Phys. 86, 487 (1999))
%   ----------------------------------------------------------------------------------------
%   In every layer j the field is a forward (+) and a backward (-) plane wave. Two 2x2 matrices
%   relate the amplitudes [E+; E-] on either side of an interface or a layer:
%
%     interface j -> k :  I_jk = 1/t_jk * [1 r_jk; r_jk 1]
%                         r_jk, t_jk = Fresnel reflection / transmission coefficients
%                         (normal incidence: r = (n_j - n_k)/(n_j + n_k), t = 2 n_j/(n_j + n_k))
%     layer j          :  L_j  = [exp(-i xi_j d_j) 0; 0 exp(+i xi_j d_j)],  xi_j = 2 pi n_j cos(theta_j) / lambda
%
%   The whole stack is the product, light side first:
%       S = I_12 L_2 I_23 L_3 ... L_(M-1) I_(M-1,M)
%   With no light coming back from the substrate, [E1+; E1-] = S [EM+; 0], so the amplitude
%   reflection coefficient is r = S21 / S11 and the power reflectance R = |S21 / S11|^2.
%
%   Oblique incidence (used for NA > 0): Snell's law fixes n_j sin(theta_j) = n_1 sin(theta_1) in
%   every layer, so cos(theta_j) = sqrt(1 - (n_1 sin(theta_1) / n_j)^2), taking the root that
%   decays into absorbing media. The Fresnel coefficients become
%       s: r = (n_j c_j - n_k c_k) / (n_j c_j + n_k c_k),   t = 2 n_j c_j / (n_j c_j + n_k c_k)
%       p: r = (n_k c_j - n_j c_k) / (n_k c_j + n_j c_k),   t = 2 n_j c_j / (n_k c_j + n_j c_k)
%   An objective illuminates the sample with a cone of angles up to asin(NA); the camera sees
%       R_NA = int (Rs + Rp)/2 w(theta) sin(theta) cos(theta) dtheta / int w(theta) sin(theta) cos(theta) dtheta
%   which shifts interference colours slightly towards the blue for large NA.
%
%   ----------------------------------------------------------------------------------------
%   Relation to TransferMatrix_packing.m / TransferMatrix_multiple.m (legacy/TransferMatrix/)
%   ----------------------------------------------------------------------------------------
%   Kept, and checked against spectra saved from the old code in tests/run_tests.m:
%     - the I_mat / L_mat matrices and R = |S21/S11|^2 of the Burkhard-Hoke-McGehee code;
%     - linear interpolation and extrapolation of n,k from the same library (tcd_materials).
%   Changed:
%     - vectorised over wavelength (one pass instead of a loop per wavelength);
%     - the materials, thicknesses and effective-medium mixtures come from a system file, not
%       from editing the function. The old code applied Maxwell-Garnett only to layers named
%       'PS-beads' and read their fill fractions with eval('vol_incl_' + i); a system file can
%       mix any two materials with any rule (see systems/README.md);
%     - oblique incidence and the NA cone average are new;
%     - the illuminant is no longer folded into R: the old output was
%       R * D65 / (12 * max(R * D65)), which mixed the lamp into the spectrum and rescaled every
%       structure differently. R is now the true reflectance, and the lamp and the eye are
%       applied once, in the colorimetry (refl_to_xyz).
%   Removed (not needed for colour, still in legacy/TransferMatrix/): the electric-field profile,
%   absorption and generation rate, the AM1.5 short-circuit current, and the incoherent thick
%   first layer (with 'Air' as the first layer, as in all our stacks, it had no effect).
%
%   See also TCD_BUILD_REFERENCE, TCD_SPECTRUM, TCD_MATERIALS.
arguments
    N {mustBeNumeric}
    d (1,:) double
    wl (1,:) double
    opts.NA (1,1) double {mustBeGreaterThanOrEqual(opts.NA, 0), mustBeLessThan(opts.NA, 1)} = 0
    opts.Angles (1,1) double {mustBeInteger, mustBePositive} = 24
    opts.Weighting (1,:) char {mustBeMember(opts.Weighting, {'uniform', 'gaussian'})} = 'uniform'
    opts.Theta double = []
    opts.Polarisation (1,:) char {mustBeMember(opts.Polarisation, {'s', 'p'})} = 's'
end
if size(N, 1) ~= numel(d)
    error('tcd:tmm', 'N has %d media (rows) but d has %d thicknesses', size(N, 1), numel(d));
end
if size(N, 2) ~= numel(wl)
    error('tcd:tmm', 'N has %d wavelengths (columns) but wl has %d', size(N, 2), numel(wl));
end
N = complex(N);

if ~isempty(opts.Theta)
    R = oblique(N, d, wl, opts.Theta, opts.Polarisation);
elseif opts.NA == 0
    R = normal(N, d, wl);
else
    % Gauss-Legendre quadrature of the cone average (same nodes as numpy's leggauss)
    tmax = asin(opts.NA);
    [x, wq] = gauss_legendre(opts.Angles);
    theta = 0.5 * tmax * (x + 1);
    w = 0.5 * tmax * wq .* sin(theta) .* cos(theta);
    if strcmp(opts.Weighting, 'gaussian'), w = w .* exp(-(theta / tmax) .^ 2); end
    R = zeros(1, numel(wl));
    for a = 1:numel(theta)
        R = R + w(a) * 0.5 * (oblique(N, d, wl, theta(a), 's') + oblique(N, d, wl, theta(a), 'p'));
    end
    R = R / sum(w);
end
end

% ---------------------------------------------------------------------------------------------
function R = normal(N, d, wl)
% Normal incidence: S = I_12 L_2 I_23 ... I_(M-1,M), kept as four 1-by-W element arrays so the
% 2x2 products run over all wavelengths at once.
[s11, s12, s21, s22] = interface(N(1, :), N(2, :), [], [], 's');
for j = 2:size(N, 1) - 1
    [s11, s12, s21, s22] = mul_layer(s11, s12, s21, s22, 2 * pi * N(j, :) * d(j) ./ wl);
    [i11, i12, i21, i22] = interface(N(j, :), N(j + 1, :), [], [], 's');
    [s11, s12, s21, s22] = mul_interface(s11, s12, s21, s22, i11, i12, i21, i22);
end
R = abs(s21 ./ s11) .^ 2;
end

function R = oblique(N, d, wl, theta0, pol)
% One angle of incidence theta0 (in the ambient) and one polarisation.
kx = N(1, :) * sin(theta0);                      % conserved tangential component n sin(theta)
C = sqrt(1 - (kx ./ N) .^ 2);                    % cos(theta_j) in every medium
flip = imag(N .* C) < 0;                         % pick the root that decays in absorbing media
C(flip) = -C(flip);
[s11, s12, s21, s22] = interface(N(1, :), N(2, :), C(1, :), C(2, :), pol);
for j = 2:size(N, 1) - 1
    [s11, s12, s21, s22] = mul_layer(s11, s12, s21, s22, 2 * pi * N(j, :) .* C(j, :) * d(j) ./ wl);
    [i11, i12, i21, i22] = interface(N(j, :), N(j + 1, :), C(j, :), C(j + 1, :), pol);
    [s11, s12, s21, s22] = mul_interface(s11, s12, s21, s22, i11, i12, i21, i22);
end
R = abs(s21 ./ s11) .^ 2;
end

function [a, b, c, d] = interface(n1, n2, c1, c2, pol)
% I = 1/t [1 r; r 1] for the interface from medium 1 into medium 2.
if isempty(c1)                                   % normal incidence
    r = (n1 - n2) ./ (n1 + n2);
    t = 2 * n1 ./ (n1 + n2);
elseif pol == 's'
    r = (n1 .* c1 - n2 .* c2) ./ (n1 .* c1 + n2 .* c2);
    t = 2 * n1 .* c1 ./ (n1 .* c1 + n2 .* c2);
else
    r = (n2 .* c1 - n1 .* c2) ./ (n2 .* c1 + n1 .* c2);
    t = 2 * n1 .* c1 ./ (n2 .* c1 + n1 .* c2);
end
a = 1 ./ t;  b = r ./ t;  c = r ./ t;  d = 1 ./ t;
end

function [s11, s12, s21, s22] = mul_layer(s11, s12, s21, s22, phase)
% S * L with L = diag(exp(-i phase), exp(+i phase)).
e1 = exp(-1i * phase);
e2 = exp(1i * phase);
s11 = s11 .* e1;  s21 = s21 .* e1;
s12 = s12 .* e2;  s22 = s22 .* e2;
end

function [t11, t12, t21, t22] = mul_interface(s11, s12, s21, s22, i11, i12, i21, i22)
% S * I, element by element over wavelength.
t11 = s11 .* i11 + s12 .* i21;  t12 = s11 .* i12 + s12 .* i22;
t21 = s21 .* i11 + s22 .* i21;  t22 = s21 .* i12 + s22 .* i22;
end

function [x, w] = gauss_legendre(n)
% Nodes and weights of n-point Gauss-Legendre quadrature on [-1, 1] (Golub-Welsch).
k = 1:n - 1;
beta = k ./ sqrt(4 * k .^ 2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)' .^ 2;
x = x';  w = w';
end
