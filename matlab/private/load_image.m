function [img, scale] = load_image(file, maxSide)
%LOAD_IMAGE RGB image as double in [0, 1], box-downsampled so its longer side <= maxSide.
[a, map] = imread(file);
if ~isempty(map), a = ind2rgb(a, map); end
if isinteger(a)
    a = double(a) / double(intmax(class(a)));
else
    a = double(a);
end
if size(a, 3) == 1, a = repmat(a, 1, 1, 3); end
a = a(:, :, 1:3);                                   % drop alpha
scale = 1;
[H, W, ~] = size(a);
if ~isempty(maxSide) && maxSide > 0 && max(H, W) > maxSide
    scale = maxSide / max(H, W);
    h = round(H * scale);  w = round(W * scale);
    Wr = box_weights(H, h);  Wc = box_weights(W, w);
    b = zeros(h, w, 3);
    for ch = 1:3
        b(:, :, ch) = Wr * a(:, :, ch) * Wc';
    end
    a = b;
end
img = a;
end

function M = box_weights(N, n)
% n-by-N area-averaging matrix: output pixel i covers input span [(i-1)N/n, iN/n).
edges = (0:n) * N / n;
rows = []; cols = []; vals = [];
for i = 1:n
    lo = edges(i);  hi = edges(i + 1);
    for r = floor(lo) + 1:ceil(hi)
        ov = min(hi, r) - max(lo, r - 1);
        if ov > 0
            rows(end + 1) = i; cols(end + 1) = r; vals(end + 1) = ov / (hi - lo); %#ok<AGROW>
        end
    end
end
M = sparse(rows, cols, vals, n, N);
end
