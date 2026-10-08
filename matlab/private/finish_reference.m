function ref = finish_reference(name, params, label, XYZ)
%FINISH_REFERENCE Fill in the colour columns of a reference from its XYZ.
%   PARAMS is a struct of columns (sweep parameters and the label column); LABEL describes the
%   label column (name, unit, title, classes, class_name, zero_name, gap).
ref.name = name;
ref.params = params;
ref.label = label;
ref.value = params.(label.name)(:);
ref.XYZ = XYZ;
ref.Jab = xyz_to_cam02ucs(XYZ);
ref.Lab = xyz_to_lab(XYZ);
c = colour_const();
ref.sRGB = min(max(srgb_codec(XYZ / 100 * c.XYZ2RGB', 'encode'), 0), 1);
ref.dE = vecnorm(ref.Jab - ref.Jab(1, :), 2, 2);
end
