"""Synthetic check of `run_tcd.py gamma`: a linear and an sRGB camera image the same field.

Run from Codes/TCD_pipeline:  python tests/check_gamma_fit.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcd import colorimetry as C  # noqa: E402

rng = np.random.default_rng(0)
scene = np.array([0.19, 0.21, 0.42])  # linear substrate intensity per channel at t = 1
exposures = [0.25, 0.5, 1.0, 1.5, 2.0]
here = Path(__file__).resolve().parents[1]
ok = True
with tempfile.TemporaryDirectory() as tmp:
    for name, encode, expect in (("linear", lambda x: x, "--gamma 1.0"),
                                 ("srgb", C.colour.models.eotf_inverse_sRGB, "--gamma srgb")):
        paths = []
        for t in exposures:
            lin = np.clip(scene * t + rng.normal(0, 0.003, (120, 160, 3)), 0, 1)
            img = (np.clip(encode(lin), 0, 1) * 255).round().astype(np.uint8)
            p = Path(tmp) / f"{name}_{t}.png"
            Image.fromarray(img).save(p)
            paths.append(str(p))
        out = subprocess.run([sys.executable, "-I", str(here / "run_tcd.py"), "gamma", "--images", *paths,
                              "--exposures", *map(str, exposures)], capture_output=True, text=True).stdout
        verdict = out.strip().splitlines()[-1]
        passed = expect in verdict
        ok &= passed
        print(f"[{'PASS' if passed else 'FAIL'}] {name} camera -> {verdict}")
sys.exit(0 if ok else 1)
