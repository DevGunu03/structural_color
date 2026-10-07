function XYZ = refl_to_xyz(R, illuminant)
%REFL_TO_XYZ Tristimulus values (perfect reflector Y = 100) of reflectance spectra.
%   R          : N-by-471 reflectances on 360:830 nm
%   illuminant : 'D65' (default), 'A', 'D50', or a CSV file of (wavelength nm, power)
%   A non-D65 illuminant is followed by CAT02 adaptation to D65, as a white-balanced camera does.
c = colour_const();
if nargin < 2, illuminant = 'D65'; end
names = c.illum.Properties.VariableNames;
if any(strcmp(names, illuminant))
    S = c.illum.(illuminant)';
else
    d = readmatrix(illuminant);
    d = d(all(isfinite(d), 2), :);
    S = interp1(d(:, 1), d(:, 2), c.wl, 'linear', 0);
end
k = 100 / sum(S .* c.cmf(:, 2)');
XYZ = k * (R .* S) * c.cmf;
if ~strcmp(illuminant, 'D65')
    white = k * S * c.cmf;
    M = c.M_CAT02;
    gain = (M * c.XYZ_D65') ./ (M * white');
    XYZ = XYZ * (M \ diag(gain) * M)';
end
end
