function g = tcd_measure_gamma(files, exposures, roi)
%TCD_MEASURE_GAMMA Camera decoding exponent from an exposure series.
%
%   g = tcd_measure_gamma({'t5.tif','t10.tif','t20.tif','t40.tif'}, [5 10 20 40])
%   g = tcd_measure_gamma(files, exposures, [x y w h])    % evenly lit region, 1-based corner
%
%   Image one field (a bare substrate is ideal) at 4-6 exposure times with auto-exposure,
%   auto-gain and auto-white-balance OFF. Pixel value v grows as t^s: s = 1 for a linear camera,
%   s ~ 0.42-0.45 for sRGB-encoded output. Returns g = 1/s per channel, the value to pass as
%   'Gamma' to TCD_MAP_IMAGE (or 'srgb' when g is about 2.2).
if nargin < 3, roi = []; end
v = nan(numel(files), 3);
for i = 1:numel(files)
    img = load_image(files{i}, 0);
    if isempty(roi)
        [H, W, ~] = size(img);
        patch = img(floor(H / 4) + 1:floor(3 * H / 4), floor(W / 4) + 1:floor(3 * W / 4), :);
    else
        patch = img(roi(2):roi(2) + roi(4) - 1, roi(1):roi(1) + roi(3) - 1, :);
    end
    px = reshape(patch, [], 3);
    clipped = mean(px < 0.02 | px > 0.985, 1) > 0.01;
    m = mean(px, 1);
    m(clipped) = nan;
    v(i, :) = m;
    fprintf('  t = %g: mean R,G,B = %s  (%s)\n', exposures(i), mat2str(m, 4), files{i});
end
s = nan(1, 3);
for ch = 1:3
    ok = isfinite(v(:, ch));
    if nnz(ok) >= 3
        p = polyfit(log(exposures(ok)), log(v(ok, ch)), 1);
        s(ch) = p(1);
    end
end
g = 1 ./ s;
gm = mean(g, 'omitnan');
fprintf('log-log slope (R,G,B) = %s  ->  decoding exponent g = %s\n', mat2str(s, 3), mat2str(g, 3));
if abs(gm - 1) < 0.25
    fprintf('=> linear camera: use ''Gamma'', 1\n');
elseif gm > 1.9 && gm < 2.6
    fprintf('=> sRGB-like encoding: use ''Gamma'', ''srgb''\n');
else
    fprintf('=> use ''Gamma'', %.2f\n', gm);
end
end
