"""Self-contained checks of the Python TCD code (no lab images needed).

Run from the repository root:  python python/tests/test_pipeline.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from tcd import colorimetry as C  # noqa: E402
from tcd.image_tcd import analyse, layer_classes  # noqa: E402
from tcd.materials import F_HEX_MONOLAYER, maxwell_garnett, nk  # noqa: E402
from tcd.reference import Reference  # noqa: E402
from tcd.tmm import reflectance, reflectance_na  # noqa: E402

FX = ROOT / "tests" / "fixtures"
fails = 0


def check(name, ok):
    global fails
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")


def fixture(name):
    d = np.loadtxt(FX / name, delimiter=",")
    return d[:, 0], d[:, 1]


# 1. transfer matrix against the original TransferMatrix_multiple/_packing.m output
wl, ref = fixture("matlab_PS_solid_film_600nm.txt")
air, ps = nk("Air", wl), nk("PS-beads", wl)
R = reflectance([air, maxwell_garnett(air, ps, 1.0), nk("SiO2-Franta", wl), nk("Si-Franta", wl)], [0, 600, 100, 0], wl)
check("TMM vs TransferMatrix_multiple.m, solid PS 600 nm", np.abs(R - ref).max() < 1e-10)
wl, ref = fixture("matlab_PS_bilayer_top50pct.txt")
layers = [air, maxwell_garnett(air, ps, F_HEX_MONOLAYER / 2), maxwell_garnett(air, ps, F_HEX_MONOLAYER),
          nk("SiO2-Franta", wl), nk("Si-Franta", wl)]
check("TMM vs TransferMatrix_packing.m, PS bilayer 50 %", np.abs(reflectance(layers, [0, 300, 300, 100, 0], wl) - ref).max() < 1e-5)
wl, ref = fixture("matlab_MoO3_100nm.txt")
R = reflectance([nk(m, wl) for m in ("Air", "aMoO3", "SiO2-Franta", "Si-Franta")], [0, 100, 100, 0], wl)
check("TMM vs TransferMatrix_multiple.m, MoO3 100 nm", np.abs(R - ref).max() < 1e-10)
check("NA -> 0 limit of the cone-averaged TMM", np.abs(reflectance_na([nk(m, wl) for m in ("Air", "aMoO3", "SiO2-Franta", "Si-Franta")],
                                                                      [0, 100, 100, 0], wl, 0.02) - R).max() < 1e-3)

# 2. CAM02-UCS against Cobeldick's sRGB_to_CAM02UCS.m (value computed with the toolbox)
jab = C.XYZ_to_cam02ucs(C.linear_to_XYZ(C.srgb_to_linear(np.array([0.34712092, 0.44539356, 0.59242291]))))
check("CAM02-UCS vs CIECAM02 toolbox", np.allclose(jab, [49.19585505, -5.70443096, -17.40028325], atol=0.01))

# 3. CIELAB of bare 100 nm SiO2/Si is physical (the old notebooks gave L* ~ 1)
moo3 = Reference.load(ROOT / "refs" / "moo3_D65.csv")
check("CIELAB L* of bare substrate", 40 < moo3.Lab[0, 0] < 60)

# 4-5. synthetic micrographs: simulated colours, unknown camera white balance, noise
gains = np.array([0.55, 0.80, 1.25])
rng = np.random.default_rng(0)


def synth(lin_colours, folder):
    img = np.zeros((80, 240, 3))
    for b, col in enumerate(lin_colours):
        img[:, b * 60:(b + 1) * 60] = col * gains
    img = np.clip(img + rng.normal(0, 0.004, img.shape), 0, 1)
    path = Path(folder) / "synthetic.png"
    Image.fromarray((img * 255).round().astype(np.uint8)).save(path)
    return path


ps_ref = Reference.load(ROOT / "refs" / "ps_D65.csv")
with tempfile.TemporaryDirectory() as tmp:
    rows = [0] + [ps_ref.find(f"{k}L") for k in (1, 2, 3)]
    res = analyse(synth(C.XYZ_to_linear(ps_ref.XYZ[rows]), tmp), ps_ref, (10, 10, 30, 30))
    cls = layer_classes(res)
    got = [int(np.bincount(cls[20:60, b * 60 + 15:(b + 1) * 60 - 15].ravel() + 1).argmax() - 1) for b in range(4)]
    check(f"synthetic PS image: layers {got}", got == [0, 1, 2, 3])

    t = [0, 150, 280, 420]
    res = analyse(synth(C.XYZ_to_linear(moo3.XYZ[t]), tmp), moo3, (10, 10, 30, 30))
    got = [float(np.nanmedian(res.value[20:60, b * 60 + 15:(b + 1) * 60 - 15])) for b in range(4)]
    check(f"synthetic MoO3 image: {got} nm", np.all(np.abs(np.array(got) - t) <= 3))

print("\nall checks passed" if fails == 0 else f"\n{fails} check(s) failed")
sys.exit(fails > 0)
