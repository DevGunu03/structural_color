"""Command-line entry point (run from the repository root).

Build a simulated reference:
    python python/run_tcd.py build-ref --system MoO3 --t-max 600 --out refs/moo3_D65.csv
    python python/run_tcd.py build-ref --system PS --max-layers 5 --out refs/ps_D65.csv

Map a micrograph (substrate region given as x,y,w,h in image pixels, or picked by clicking):
    python python/run_tcd.py map --ref refs/ps_D65.csv --image my_image.png --interactive --out results/my_image
    python python/run_tcd.py map --ref refs/ps_D65.csv --image my_image.png --substrate 186,433,40,40 --out results/my_image
    python python/run_tcd.py map --ref refs/ps_D65.csv --image my_image.png --substrate-colour 0.157,0.263,0.459 \
                                 --out results/my_image
    python python/run_tcd.py map --ref refs/moo3_D65.csv --image flakes.png --substrate auto --out results/flakes

Measure the camera's decoding exponent (same field, auto-exposure/gain/white balance off):
    python python/run_tcd.py gamma --images t5.tif t10.tif t20.tif t40.tif --exposures 5 10 20 40
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


def _floats(s: str) -> list[float]:
    return [float(v) for v in s.split(",")]


def cmd_build_ref(args):
    import matplotlib
    matplotlib.use("Agg")
    from tcd.plots import reference_figure
    from tcd.reference import build_moo3, build_ps

    if args.system.lower() == "moo3":
        ref = build_moo3(t_max=args.t_max, step=args.step, oxide_nm=args.oxide, moo3=args.material or "aMoO3",
                         illuminant=args.illuminant, na=args.na)
    else:
        ref = build_ps(max_layers=args.max_layers, packing_step=args.packing_step, bead_nm=args.bead,
                       oxide_nm=args.oxide, model=args.model, ps=args.material or "PS-beads",
                       illuminant=args.illuminant, na=args.na)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ref.save(out, with_spectra=args.spectra)
    reference_figure(ref, out.with_suffix(".png"))
    print(f"{ref.system} reference: {len(ref)} rows -> {out} (+ .json, .png)")


def _resolve_roi(spec: str | None, img, scale, ref, label, interactive, size=40):
    from tcd.roi import dominant_colour, select_roi_interactive, suggest_roi

    if spec is None:
        if interactive:
            return select_roi_interactive(img, f"Select a {label} region", scale)
        raise SystemExit(f"Give a region for '{label}' (x,y,w,h | auto | colour:r,g,b) or use --interactive")
    if spec == "auto":
        return suggest_roi(img, dominant_colour(img), size=size, scale=scale)
    if spec.startswith("colour:") or spec.startswith("color:"):
        return suggest_roi(img, _floats(spec.split(":", 1)[1]), size=size, scale=scale)
    x, y, w, h = (int(v) for v in _floats(spec))
    return (x, y, w, h)


def cmd_map(args):
    import matplotlib
    if not args.interactive:
        matplotlib.use("Agg")
    from tcd.image_tcd import analyse, layer_classes, load_image, summary_rows
    from tcd.plots import result_figure
    from tcd.reference import Reference

    ref = Reference.load(args.ref)
    if ref.system == "MoO3" and args.t_max:
        ref = ref.subset(ref.value <= args.t_max)
    if ref.system == "PS" and args.max_layers:
        ref = ref.subset(ref.params["layers"] <= args.max_layers)

    preview, scale = load_image(args.image, args.max_side)
    sub_spec = args.substrate
    if args.substrate_colour:
        sub_spec = "colour:" + args.substrate_colour
    substrate = _resolve_roi(sub_spec, preview, scale, ref, "bare-substrate", args.interactive, args.roi_size)
    anchors = []
    for a in args.anchor or []:
        label, spec = a.split("@", 1)
        roi = _resolve_roi(spec, preview, scale, ref, label, args.interactive, args.anchor_size)
        anchors.append((roi, ref.find(label), label))

    gamma = args.gamma if args.gamma == "srgb" else float(args.gamma)
    res = analyse(args.image, ref, substrate, anchors, max_side=args.max_side, smooth_sigma=args.sigma,
                  max_residual=args.max_residual, correction_mode=args.correction, gamma=gamma)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result_figure(res, out.with_suffix(".png"), title=f"{Path(args.image).name}  |  {ref.meta.get('stack', '')}")
    rows = summary_rows(res)
    with open(out.with_name(out.name + "_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    maps = dict(tcd=res.tcd, value=res.value, residual=res.residual)
    if ref.system == "PS":
        maps["layers"] = layer_classes(res)
    if res.alt_value is not None:
        maps.update(alt_value=res.alt_value, alt_residual=res.alt_residual)
    np.savez_compressed(out.with_name(out.name + "_maps.npz"), **maps)
    info = dict(image=str(args.image), reference=str(args.ref), anchors=[a.__dict__ for a in res.anchors],
                correction_matrix=res.correction.tolist(), settings=res.settings)
    out.with_name(out.name + "_run.json").write_text(json.dumps(info, indent=2, default=str))

    print(f"{Path(args.image).name}: substrate ROI {substrate}, anchors {[(a.name, a.roi) for a in res.anchors[1:]]}")
    print(f"  median residual (assigned) = {np.nanmedian(res.residual[res.assigned]):.2f}, "
          f"unassigned = {100 * (~res.assigned).mean():.1f} %")
    for r in rows:
        if r["fraction"] > 0.005:
            print(f"  {r['item']:>12s}: {100 * r['fraction']:5.1f} %   median residual {r['median_residual']:.2f}")
    print(f"  -> {out.with_suffix('.png')}")


def cmd_gamma(args):
    from tcd.camera import fit_decoding_gamma, region_means

    if len(args.images) != len(args.exposures):
        raise SystemExit("--images and --exposures need the same number of entries")
    roi = tuple(int(v) for v in _floats(args.roi)) if args.roi else None
    values = region_means(args.images, roi)
    for p, t, v in zip(args.images, args.exposures, values):
        print(f"  t = {t:g}: mean R,G,B = {np.round(v, 4)}  ({Path(p).name})")
    slopes, g = fit_decoding_gamma(values, np.array(args.exposures))
    print(f"log-log slope (R,G,B) = {np.round(slopes, 3)}  ->  decoding exponent g = {np.round(g, 2)}")
    g_mean = float(np.nanmean(g))
    verdict = "linear camera: use --gamma 1.0" if abs(g_mean - 1) < 0.25 else \
        "sRGB-like encoding: use --gamma srgb" if 1.9 < g_mean < 2.6 else f"use --gamma {g_mean:.2f}"
    print(f"=> {verdict}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-ref", help="simulate a colour reference")
    b.add_argument("--system", required=True, choices=["MoO3", "PS", "moo3", "ps"])
    b.add_argument("--out", required=True)
    b.add_argument("--illuminant", default="D65", help="D65, A, D50, ... or a CSV of a measured lamp spectrum")
    b.add_argument("--na", type=float, default=0.0, help="objective NA for cone-averaged reflectance (0 = normal)")
    b.add_argument("--oxide", type=float, default=100.0, help="SiO2 thickness (nm)")
    b.add_argument("--material", help="override MoO3 / PS optical constants (library name or n,k file)")
    b.add_argument("--t-max", type=float, default=600.0)
    b.add_argument("--step", type=float, default=1.0)
    b.add_argument("--max-layers", type=int, default=5)
    b.add_argument("--packing-step", type=float, default=0.01)
    b.add_argument("--bead", type=float, default=300.0, help="PS bead diameter (nm)")
    b.add_argument("--model", default="slab", choices=["slab", "hcp"])
    b.add_argument("--spectra", action="store_true", help="also save the reflectance spectra")
    b.set_defaults(func=cmd_build_ref)

    m = sub.add_parser("map", help="thickness / layer map of a micrograph")
    m.add_argument("--ref", required=True)
    m.add_argument("--image", required=True)
    m.add_argument("--out", required=True, help="output path stem")
    m.add_argument("--substrate", help="x,y,w,h | auto | colour:r,g,b")
    m.add_argument("--substrate-colour", help="r,g,b (0-1) of the bare substrate; region is found automatically")
    m.add_argument("--anchor", action="append",
                   help="extra calibration region of known structure, e.g. 1L@120,40,30,30 or 1L@colour:r,g,b")
    m.add_argument("--interactive", action="store_true", help="click regions that were not given")
    m.add_argument("--gamma", default="1.0",
                   help="camera decoding: 'srgb' or a number g (linear = g**value); default 1.0 = linear camera")
    m.add_argument("--roi-size", type=int, default=40, help="side (px) of auto-placed substrate region")
    m.add_argument("--anchor-size", type=int, default=12, help="side (px) of auto-placed anchor regions")
    m.add_argument("--correction", default="auto", choices=["auto", "diagonal", "matrix"])
    m.add_argument("--max-residual", type=float, default=15.0)
    m.add_argument("--sigma", type=float, default=1.2, help="Gaussian smoothing (px)")
    m.add_argument("--max-side", type=int, default=1600, help="downsample longer side to this (0 = off)")
    m.add_argument("--t-max", type=float, help="MoO3: restrict reference to thickness <= t-max (prior)")
    m.add_argument("--max-layers", type=int, help="PS: restrict reference to <= this many layers")
    m.set_defaults(func=cmd_map)

    g = sub.add_parser("gamma", help="measure the camera's decoding exponent from an exposure series")
    g.add_argument("--images", nargs="+", required=True, help="same field, only exposure time changed")
    g.add_argument("--exposures", nargs="+", type=float, required=True, help="exposure times, same order")
    g.add_argument("--roi", help="x,y,w,h of an evenly lit region (default: central half of the frame)")
    g.set_defaults(func=cmd_gamma)

    args = p.parse_args(argv)
    if getattr(args, "max_side", None) == 0:
        args.max_side = None
    args.func(args)


if __name__ == "__main__":
    main()
