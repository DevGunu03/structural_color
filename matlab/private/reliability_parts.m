function [fit, cons] = reliability_parts(value, residual, ref, tol, sigma)
%RELIABILITY_PARTS Fit and neighbour consistency of a thickness / layer map (see tcd_map_image).
%   fit  = P(colour error >= residual) for a 3-D Gaussian error with SIGMA per axis
%        = erfc(x / sqrt 2) + sqrt(2 / pi) x exp(-x^2 / 2), x = residual / sigma
%   cons = share of assigned pixels in the 5 x 5 window that agree with the centre pixel
%   Same as chi3_sf and consistency in python/tcd/reliability.py.
x = residual / sigma;
fit = erfc(x / sqrt(2)) + sqrt(2 / pi) * x .* exp(-x .^ 2 / 2);
fit(isnan(value)) = nan;
[H, W] = size(value);
r = 2;
padded = nan(H + 2 * r, W + 2 * r);
padded(r + 1:r + H, r + 1:r + W) = value;
agree = zeros(H, W);
count = zeros(H, W);
for dy = 0:2 * r
    for dx = 0:2 * r
        nb = padded(dy + 1:dy + H, dx + 1:dx + W);
        ok = ~isnan(nb);
        count = count + ok;
        agree = agree + (ok & same_structure(ref, nb, value, tol));
    end
end
cons = agree ./ max(count, 1);
cons(isnan(value)) = nan;
end
