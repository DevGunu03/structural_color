"""Colour-based thickness / layer-number mapping (TCD) of thin films from optical micrographs.

system.py   system files (systems/*.jsonc) -> candidate structures and their layer stacks
tmm.py      transfer-matrix reflectance (normal incidence and objective NA)
materials.py  optical constants: n,k library, files, Cauchy, Sellmeier, effective media
colorimetry.py  spectra -> XYZ -> CAM02-UCS / CIELAB / sRGB
reference.py  simulated colour references (build, save, load)
image_tcd.py  micrograph -> calibrated colours -> thickness / layer maps
"""
