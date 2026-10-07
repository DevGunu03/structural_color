function out = tcd_colour(conversion, in, illuminant)
%TCD_COLOUR Colour conversions used by the TCD code (rows = colours).
%
%   Jab = tcd_colour('srgb2cam02ucs', [0.35 0.45 0.59])   % encoded sRGB (0-1) -> CAM02-UCS J'a'b'
%   Jab = tcd_colour('xyz2cam02ucs', XYZ)                 % XYZ on the 0-100 scale, D65
%   Lab = tcd_colour('xyz2lab', XYZ)                      % CIELAB, D65 white
%   rgb = tcd_colour('xyz2srgb', XYZ)                     % encoded sRGB, clipped to [0, 1]
%   XYZ = tcd_colour('refl2xyz', R, 'D65')                % reflectance spectra (N-by-471, 360:830 nm)
%
%   CAM02-UCS uses the viewing conditions of Cobeldick's sRGB_to_CAM02UCS.m
%   (L_A = 64/(5 pi) cd/m^2, Y_b = 20, average surround), so 'srgb2cam02ucs' reproduces the
%   earlier step 6_CIEUCS_from_sRGB.m without the CIECAM02 toolbox.
c = colour_const();
switch lower(conversion)
    case 'srgb2cam02ucs'
        out = xyz_to_cam02ucs(100 * srgb_codec(in, 'decode') * c.RGB2XYZ');
    case 'xyz2cam02ucs'
        out = xyz_to_cam02ucs(in);
    case 'xyz2lab'
        out = xyz_to_lab(in);
    case 'xyz2srgb'
        out = min(max(srgb_codec(in / 100 * c.XYZ2RGB', 'encode'), 0), 1);
    case 'refl2xyz'
        if nargin < 3, illuminant = 'D65'; end
        out = refl_to_xyz(in, illuminant);
    otherwise
        error('tcd:colour', 'Unknown conversion "%s"', conversion);
end
end
