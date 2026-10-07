function ref = finish_reference(system, params, value, XYZ)
%FINISH_REFERENCE Fill in the colour columns of a reference from its XYZ.
ref.system = system;
ref.params = params;
ref.value = value(:);
ref.XYZ = XYZ;
ref.Jab = xyz_to_cam02ucs(XYZ);
ref.Lab = xyz_to_lab(XYZ);
c = colour_const();
ref.sRGB = min(max(srgb_codec(XYZ / 100 * c.XYZ2RGB', 'encode'), 0), 1);
ref.dE = vecnorm(ref.Jab - ref.Jab(1, :), 2, 2);
end
