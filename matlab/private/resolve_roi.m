function roi = resolve_roi(spec, img, scale, label, defaultSize)
%RESOLVE_ROI Region [x y w h] in original-image pixels (x, y = 1-based top-left corner).
%   spec : [x y w h] | 'click' | 'auto' (dominant colour) | [r g b] colour in 0-1
if isnumeric(spec) && numel(spec) == 4
    roi = double(spec(:)');
elseif isnumeric(spec) && numel(spec) == 3
    roi = suggest(img, spec(:)', defaultSize, scale);
elseif strcmpi(spec, 'auto')
    roi = suggest(img, dominant_colour(img), defaultSize, scale);
elseif strcmpi(spec, 'click')
    f = figure('Name', ['Select ' label]);
    image(img); axis image off;
    title(sprintf('Click two opposite corners of a %s region', label));
    [x, y] = ginput(2);
    close(f);
    if numel(x) < 2, error('tcd:roi', 'Region selection cancelled'); end
    roi = [floor(min(x) / scale) + 1, floor(min(y) / scale) + 1, ...
           max(round(abs(diff(x)) / scale), 1), max(round(abs(diff(y)) / scale), 1)];
else
    error('tcd:roi', 'Region for %s must be [x y w h], [r g b], ''auto'' or ''click''', label);
end
end

function roi = suggest(img, target, sz, scale)
% Square window whose pixels are, on average, closest to TARGET (sRGB 0-1).
d = sqrt(sum((img - reshape(target, 1, 1, 3)) .^ 2, 3));
s = max(3, round(sz * scale));
B = conv2(d, ones(s) / s ^ 2, 'valid');
[~, i] = min(B(:));
[r, c] = ind2sub(size(B), i);
roi = [floor((c - 1) / scale) + 1, floor((r - 1) / scale) + 1, sz, sz];
end

function col = dominant_colour(img)
% Most common colour (16-level 3-D histogram mode): the substrate in sparse-flake images.
bins = 16;
px = reshape(img, [], 3);
q = min(floor(px * bins), bins - 1);
code = q(:, 1) * bins ^ 2 + q(:, 2) * bins + q(:, 3);
m = mode(code);
col = mean(px(code == m, :), 1);
end
