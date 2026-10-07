"""Map the bundled example micrographs with the Python code.

Run from the repository root:  python examples/run_examples.py
Results go to results/examples/ (figures, summaries, per-pixel maps).
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMG = ROOT / "examples" / "images"
OUT = ROOT / "results" / "examples"

# Substrate colours (sRGB, 0-1) of the two PS microscopes, read off bare-substrate areas;
# the MoO3 images use the dominant background colour instead.
CASES = [
    ("ps_chennai_metco", "ps_D65.csv", ["--substrate-colour", "0.157,0.263,0.459"]),
    ("ps_olympus", "ps_D65.csv", ["--substrate-colour", "0.580,0.561,0.333"]),
    ("moo3_region_1", "moo3_D65.csv", ["--substrate", "auto"]),
    ("moo3_region_3", "moo3_D65.csv", ["--substrate", "auto"]),
    ("moo3_50x2", "moo3_D65.csv", ["--substrate", "auto"]),
]

for name, ref, extra in CASES:
    cmd = [sys.executable, str(ROOT / "python" / "run_tcd.py"), "map", "--ref", str(ROOT / "refs" / ref),
           "--image", str(IMG / f"{name}.png"), "--out", str(OUT / name), *extra]
    subprocess.run(cmd, check=True)
