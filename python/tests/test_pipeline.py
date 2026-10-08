"""Self-contained checks of the Python TCD code (no lab images needed).

Run from the repository root:  python python/tests/test_pipeline.py
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from tcd import colorimetry as C  # noqa: E402
from tcd.expr import ExpressionError, evaluate, names_in  # noqa: E402
from tcd.image_tcd import analyse, layer_classes  # noqa: E402
from tcd.materials import F_HEX_MONOLAYER, bruggeman, maxwell_garnett, nk, sellmeier  # noqa: E402
from tcd.reference import Reference, build_reference  # noqa: E402
from tcd.system import SYSTEMS_DIR, SystemError, load_system  # noqa: E402
from tcd.tmm import reflectance, reflectance_na  # noqa: E402

FX = ROOT / "tests" / "fixtures"
WL = C.WAVELENGTHS
fails = 0


def check(name, ok):
    global fails
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")


def raises(fn, exc=Exception, text=""):
    try:
        fn()
    except exc as e:
        return text.lower() in str(e).lower()
    return False


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
moo3_stack = [nk(m, wl) for m in ("Air", "aMoO3", "SiO2-Franta", "Si-Franta")]
R = reflectance(moo3_stack, [0, 100, 100, 0], wl)
check("TMM vs TransferMatrix_multiple.m, MoO3 100 nm", np.abs(R - ref).max() < 1e-10)
check("NA -> 0 limit of the cone-averaged TMM", np.abs(reflectance_na(moo3_stack, [0, 100, 100, 0], wl, 0.02) - R).max() < 1e-3)
check("zero-thickness layer changes nothing", np.abs(reflectance(moo3_stack, [0, 0, 100, 0], wl) - reflectance(
    [moo3_stack[0], moo3_stack[2], moo3_stack[3]], [0, 100, 0], wl)).max() < 1e-12)

# 2. colour conversions
jab = C.XYZ_to_cam02ucs(C.linear_to_XYZ(C.srgb_to_linear(np.array([0.34712092, 0.44539356, 0.59242291]))))
check("CAM02-UCS vs CIECAM02 toolbox", np.allclose(jab, [49.19585505, -5.70443096, -17.40028325], atol=0.01))
moo3 = Reference.load(ROOT / "refs" / "moo3_D65.csv")
check("CIELAB L* of bare substrate", 40 < moo3.Lab[0, 0] < 60)

# 3. expression language (same rules as matlab/private/expr_eval.m)
cases = {"1 + 2 * 3": 7, "-2^2": -4, "2^-1": 0.5, "2**3": 8, "(1 + 2) * 3": 9, "max(1, min(5, 3))": 3,
         "round(2.5) + round(-2.5)": 1, "3 >= 2": 1, "2 == 3": 0, "sqrt(16) / 2e0": 2, "pi": np.pi,
         "1 - 2 - 3": -4, "8 / 4 / 2": 1, "2^3^2": 512, ".5 + 1e-3": 0.501}
got = {e: float(evaluate(e)) for e in cases}
check("expressions: precedence, functions, comparisons", all(np.isclose(got[e], v) for e, v in cases.items()))
check("expressions: arrays broadcast", np.allclose(evaluate("2 * t + c", {"t": np.arange(3.0), "c": 1}), [1, 3, 5]))
check("expressions: unknown names and code are refused",
      raises(lambda: evaluate("t + 1"), ExpressionError, "unknown name")
      and raises(lambda: evaluate("__import__('os')"), ExpressionError)
      and raises(lambda: evaluate("2 +"), ExpressionError) and raises(lambda: evaluate("max(1)"), ExpressionError))
check("expressions: names_in", names_in("max(layers - 1, 0) * f_hex + pi") == {"layers", "f_hex"})

# 4. material models
a, b = nk("Air", WL), nk("PS-beads", WL)
check("Maxwell-Garnett and Bruggeman reduce to host / inclusion at f = 0 / 1",
      max(np.abs(maxwell_garnett(a, b, 0) - a).max(), np.abs(bruggeman(a, b, 1) - b).max(),
          np.abs(bruggeman(a, b, 0) - a).max()) < 1e-12)
check("Bruggeman is symmetric in host and inclusion", np.abs(bruggeman(a, b, 0.3) - bruggeman(b, a, 0.7)).max() < 1e-12)
n_d = sellmeier(np.array([587.6]), [0.6961663, 0.4079426, 0.8974794], [0.0684043**2, 0.1162414**2, 9.896161**2])
check("Sellmeier: fused silica n_d = 1.4585", abs(n_d[0].real - 1.4585) < 2e-4)

# 5. system files
systems = sorted(SYSTEMS_DIR.glob("*.jsonc"))
loaded = {}
for f in systems:
    try:
        loaded[f.stem] = load_system(f)
    except Exception as e:  # noqa: BLE001
        print("   ", f.name, e)
check(f"all {len(systems)} bundled system files load", len(loaded) == len(systems))
for name in ("moo3", "ps", "ps_hcp", "sio2_on_si", "graphene", "bragg_mirror"):
    committed = Reference.load(ROOT / "refs" / f"{name}_D65.csv")
    rebuilt = build_reference(loaded[name])
    same_rows = all(np.allclose(rebuilt.params[k], committed.params[k], atol=1e-6) for k in committed.params)
    rel = np.max(np.abs(rebuilt.XYZ - committed.XYZ) / np.maximum(committed.XYZ, 1e-3))
    check(f"systems/{name}.jsonc reproduces refs/{name}_D65.csv (max rel XYZ diff {rel:.1e})", same_rows and rel < 2e-5)
s = load_system("ps", {"oxide_nm": "285", "max_layers": "2"})
check("--set overrides constants", s.constants["oxide_nm"] == 285 and s.rows["layers"].max() == 2)
two = load_system("bragg_mirror", {"pairs": "2"}).at([{"lambda0": 550}])
N, d, names = two.stack(0)
explicit = reflectance(N, d, WL)
manual = [nk("Air", WL)] + [nk("TiO2-Franta", WL), nk("SiO2-Franta", WL)] * 2 + [nk("Si-Franta", WL)]
check("a repeated group equals the layers written out", len(N) == 6 and np.abs(
    explicit - reflectance(manual, [0] + [550 / 4 / 2.35, 550 / 4 / 1.46] * 2 + [0], WL)).max() < 1e-12)
base = json.loads(json.dumps(loaded["moo3"].definition))
bad = [("unknown key", dict(base, thicknes=1), "unknown key"),
       ("unknown material", dict(base, substrate="Unobtainium"), "not in the n,k library"),
       ("fraction > 1", dict(base, materials={"m": {"ema": "bruggeman", "host": "Air", "inclusion": "aMoO3",
                                                    "fraction": "thickness_nm / 100"}},
                             layers=[{"material": "m", "thickness": 5}]), "within [0, 1]"),
       ("bad expression", dict(base, layers=[{"material": "aMoO3", "thickness": "thickness_nm +"}]), "thickness"),
       ("unknown --set", None, "not a constant")]
ok = all(raises(lambda d=d: load_system(d), SystemError, t) for _, d, t in bad[:-1])
ok &= raises(lambda: load_system("moo3", {"oxide": "285"}), SystemError, "not a constant")
check("system file mistakes give clear errors", ok)

# 6. physics sanity: monolayer graphene contrast peaks at the oxide's reflectance minimum
g = load_system("graphene", {"oxide_nm": "300"}).at([{"layers": 0}, {"layers": 1}])
R0, R1 = g.spectra()
vis = (WL >= 450) & (WL <= 700)
contrast = (R0 - R1) / R0
peak_wl, rmin_wl = WL[vis][np.argmax(contrast[vis])], WL[vis][np.argmin(R0[vis])]
check(f"graphene on 300 nm SiO2: contrast {contrast[vis].max():.3f} at {peak_wl:.0f} nm (R min at {rmin_wl:.0f} nm)",
      0.08 < contrast[vis].max() < 0.16 and abs(peak_wl - rmin_wl) < 15)

# 7. reference files: round trip, and references written before system files
with tempfile.TemporaryDirectory() as tmp:
    gref = build_reference(loaded["graphene"])
    gref.save(Path(tmp) / "g.csv")
    back = Reference.load(Path(tmp) / "g.csv")
    check("reference save/load keeps the label description", back.label == gref.label and back.find("3L") == 3)
    meta = json.loads((ROOT / "refs" / "ps_D65.json").read_text())
    legacy = {"stack": meta["stack"], "illuminant": "D65", "NA": 0, "system": "PS", "columns": meta["columns"]}
    (Path(tmp) / "old.json").write_text(json.dumps(legacy))
    (Path(tmp) / "old.csv").write_text((ROOT / "refs" / "ps_D65.csv").read_text())
    old = Reference.load(Path(tmp) / "old.csv")
    check("old-style PS reference still loads", old.classes and old.label["name"] == "eff_layers")

# 8. synthetic micrographs: simulated colours, unknown camera white balance, noise
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


def band_modes(cls):
    return [int(np.bincount(cls[20:60, b * 60 + 15:(b + 1) * 60 - 15].ravel() + 1).argmax() - 1) for b in range(4)]


ps_ref = Reference.load(ROOT / "refs" / "ps_D65.csv")
graphene = Reference.load(ROOT / "refs" / "graphene_D65.csv")
with tempfile.TemporaryDirectory() as tmp:
    rows = [0] + [ps_ref.find(f"{k}L") for k in (1, 2, 3)]
    res = analyse(synth(C.XYZ_to_linear(ps_ref.XYZ[rows]), tmp), ps_ref, (10, 10, 30, 30))
    got = band_modes(layer_classes(res))
    check(f"synthetic PS image: layers {got}", got == [0, 1, 2, 3])

    t = [0, 150, 280, 420]
    res = analyse(synth(C.XYZ_to_linear(moo3.XYZ[t]), tmp), moo3, (10, 10, 30, 30))
    got = [float(np.nanmedian(res.value[20:60, b * 60 + 15:(b + 1) * 60 - 15])) for b in range(4)]
    check(f"synthetic MoO3 image: {got} nm", np.all(np.abs(np.array(got) - t) <= 3))

    truth = [0, 1, 3, 6]
    res = analyse(synth(C.XYZ_to_linear(graphene.XYZ[truth]), tmp), graphene, (10, 10, 30, 30))
    got = band_modes(layer_classes(res))
    check(f"synthetic graphene image: layers {got}", got == truth)

print("\nall checks passed" if fails == 0 else f"\n{fails} check(s) failed")
sys.exit(fails > 0)
