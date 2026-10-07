function out = gauss_smooth(img, sigma)
%GAUSS_SMOOTH Separable Gaussian blur per channel, reflect padding, kernel radius round(4*sigma)
%   (the same as scipy.ndimage.gaussian_filter used by the Python code).
if sigma <= 0, out = img; return; end
r = floor(4 * sigma + 0.5);
x = -r:r;
k = exp(-x .^ 2 / (2 * sigma ^ 2));
k = k / sum(k);
[H, W, C] = size(img);
ir = reflect_index(H, r);
ic = reflect_index(W, r);
out = zeros(H, W, C);
for ch = 1:C
    P = img(ir, ic, ch);
    out(:, :, ch) = conv2(k', k, P, 'valid');
end
end

function idx = reflect_index(n, r)
core = 1:n;
lo = r:-1:1;
hi = n:-1:n - r + 1;
idx = [lo, core, hi];
idx = min(max(idx, 1), n);         % images narrower than the kernel
end
