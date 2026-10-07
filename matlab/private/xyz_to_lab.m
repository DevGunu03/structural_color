function Lab = xyz_to_lab(XYZ)
%XYZ_TO_LAB CIE 1976 L*a*b* of N-by-3 XYZ (0-100 scale), D65 white.
c = colour_const();
r = XYZ ./ c.XYZ_D65;
e = 216 / 24389;  kap = 24389 / 27;
f = r .^ (1 / 3);
lo = r <= e;
f(lo) = (kap * r(lo) + 16) / 116;
Lab = [116 * f(:, 2) - 16, 500 * (f(:, 1) - f(:, 2)), 200 * (f(:, 2) - f(:, 3))];
end
