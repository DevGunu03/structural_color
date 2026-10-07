function y = srgb_codec(x, direction)
%SRGB_CODEC sRGB transfer function: 'decode' (encoded 0-1 -> linear) or 'encode'.
switch direction
    case 'decode'
        y = x / 12.92;
        hi = x > 0.04045;
        y(hi) = ((x(hi) + 0.055) / 1.055) .^ 2.4;
    case 'encode'
        y = 12.92 * x;
        hi = x > 0.0031308;
        y(hi) = 1.055 * x(hi) .^ (1 / 2.4) - 0.055;
end
end
