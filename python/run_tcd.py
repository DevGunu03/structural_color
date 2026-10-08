"""Command-line entry point (run from the repository root).

The workflow for any sample:
  1. describe it in a system file          systems/my_sample.jsonc (see systems/README.md)
  2. check single structures               python python/run_tcd.py spectrum --system systems/my_sample.jsonc --at t=100
  3. build the colour reference            python python/run_tcd.py build-ref --system systems/my_sample.jsonc
  4. map a micrograph with it              python python/run_tcd.py map --ref refs/my_sample_D65.csv --image img.png --interactive

Examples with the bundled systems:
    python python/run_tcd.py build-ref --system moo3                       # -> refs/moo3_D65.csv
    python python/run_tcd.py build-ref --system ps --set oxide_nm=285      # -> refs/ps_oxide_nm285_D65.csv
    python python/run_tcd.py spectrum --system ps --at layers=1,packing=1 --at layers=2,packing=0.5
    python python/run_tcd.py map --ref refs/ps_D65.csv --image my_image.png --substrate 186,433,40,40 --out results/my_image
    python python/run_tcd.py map --ref refs/moo3_D65.csv --image flakes.png --substrate auto --out results/flakes
    python python/run_tcd.py materials                                     # n,k library and bundled systems

Measure the camera's decoding exponent (same field, auto-exposure/gain/white balance off):
    python python/run_tcd.py gamma --images t5.tif t10.tif t20.tif t40.tif --exposures 5 10 20 40
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path(__file__).resolve().parents[1]


def _floats(s: str) -> list[float]:
    return [float(v) for v in s.split(",")]


def _assignments(items: list[str]) -> list[dict[str, float]]:
    """['layers=2,packing=0.5'] -> [{'layers': 2.0, 'packing': 0.5}]."""
    out = []
    for item in items:
        d = {}
        for part in item.split(","):
            if "=" not in part:
                raise SystemExit(f"--at expects name=value[,name=value...], got '{item}'")
            k, v = part.split("=", 1)
            d[k.strip()] = float(v)
        out.append(d)
    return out


def _default_ref_path(system, overrides: dict, illuminant: str, na: float) -> Path:
    tag = "".join(f"_{k}{v}" for k, v in overrides.items())
    tag += f"_NA{na:g}" if na else ""
    ill = Path(illuminant).stem if Path(illuminant).suffix else illuminant
    return ROOT / "refs" / re.sub(r"[^\w.-]", "_", f"{system.name}{tag}_{ill}.csv")


def cmd_build_ref(args):
    import matplotlib
    matplotlib.use("Agg")
    from tcd.plots import reference_figure
    from tcd.reference import ambiguity, build_reference
    from tcd.system import load_system, parse_overrides

    overrides = parse_overrides(args.set)
    system = load_system(args.system, overrides)
    illuminant = args.illuminant or system.optics["illuminant"]
    na = system.optics["na"] if args.na is None else args.na
    print(f"{system.name}: {system.n_rows} candidate structures")
    print(f"  stack: {system.describe()}")
    ref = build_reference(system, illuminant, na)
    out = Path(args.out) if args.out else _default_ref_path(system, overrides, illuminant, na)
    ref.save(out, with_spectra=args.spectra)
    reference_figure(ref, out.with_suffix(".png"))
    amb = ambiguity(ref)
    v = ref.value
    print(f"  label '{ref.label['name']}' from {np.nanmin(v):g} to {np.nanmax(v):g} {ref.label['unit']}; "
          f"{100 * np.mean(amb < 2):.0f} % of the candidates have a look-alike (< 2 dE) more than "
          f"{ref.label['gap']:g} {ref.label['unit']} away")
    print(f"  -> {out} (+ .json, reference sheet .png)")


def cmd_spectrum(args):
    import matplotlib
    matplotlib.use("Agg")
    from tcd import colorimetry as C
    from tcd.plots import spectra_figure
    from tcd.system import load_system, parse_overrides

    overrides = parse_overrides(args.set)
    system = load_system(args.system, overrides)
    if args.at:
        sub = system.at(_assignments(args.at))
    else:
        rows = sorted({0, system.n_rows // 2, system.n_rows - 1})
        sub = system.at([{k: float(system.rows[k][i]) for k in system.parameters
                          if not isinstance(system.grids[0][k], str)} for i in rows])
    na = system.optics["na"] if args.na is None else args.na
    illuminant = args.illuminant or system.optics["illuminant"]
    R = sub.spectra(na=na)
    XYZ = C.reflectance_to_XYZ(R, illuminant)
    Jab, rgb = C.XYZ_to_cam02ucs(XYZ), C.XYZ_to_srgb(XYZ)
    names = []
    for i in range(sub.n_rows):
        params = ", ".join(f"{k}={sub.rows[k][i]:g}" for k in sub.rows)
        names.append(params)
        N, d, layer_names = sub.stack(i)
        print(f"{params}")
        print("    stack: " + " | ".join(f"{n} {t:g} nm" if t else n for n, t in zip(layer_names, d)))
        print(f"    XYZ = {np.round(XYZ[i], 3)}, J'a'b' = {np.round(Jab[i], 2)}, sRGB = {np.round(rgb[i], 3)}, "
              f"dE from first = {np.linalg.norm(Jab[i] - Jab[0]):.2f}")
    out = Path(args.out) if args.out else ROOT / "results" / f"{system.name}_spectra.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    spectra_figure(R, rgb, names, f"{system.name}: {sub.describe()}", out)
    if args.csv:
        header = "wavelength_nm," + ",".join(f'"{n}"' for n in names)
        np.savetxt(out.with_suffix(".csv"), np.column_stack([C.WAVELENGTHS, R.T]), delimiter=",", fmt="%.6g",
                   header=header, comments="")
    print(f"-> {out}")


def _resolve_roi(spec: str | None, img, scale, label, interactive, size=40):
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


def _pick_image() -> str:
    try:
        import tkinter
        from tkinter import filedialog
    except ImportError:
        raise SystemExit("Give --image (no file dialog available)") from None
    root = tkinter.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(title="Select a micrograph",
                                      filetypes=[("Images", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp"), ("All", "*")])
    root.destroy()
    if not path:
        raise SystemExit("No image selected")
    return path


def cmd_map(args):
    import matplotlib
    if not args.interactive:
        matplotlib.use("Agg")
    from tcd.image_tcd import analyse, layer_classes, load_image, summary_rows
    from tcd.plots import result_figure
    from tcd.reference import Reference

    ref = Reference.load(args.ref)
    if args.max_value is not None:
        ref = ref.subset(ref.value <= args.max_value)
    image = args.image or _pick_image()
    out = Path(args.out) if args.out else ROOT / "results" / Path(image).stem

    preview, scale = load_image(image, args.max_side)
    sub_spec = "colour:" + args.substrate_colour if args.substrate_colour else args.substrate
    substrate = _resolve_roi(sub_spec, preview, scale, "bare-substrate", args.interactive, args.roi_size)
    anchors = []
    for a in args.anchor or []:
        label, spec = a.split("@", 1)
        roi = _resolve_roi(spec, preview, scale, label, args.interactive, args.anchor_size)
        anchors.append((roi, ref.find(label), label))

    gamma = args.gamma if args.gamma == "srgb" else float(args.gamma)
    res = analyse(image, ref, substrate, anchors, max_side=args.max_side, smooth_sigma=args.sigma,
                  max_residual=args.max_residual, correction_mode=args.correction, gamma=gamma,
                  ambiguity_gap=args.gap)

    out.parent.mkdir(parents=True, exist_ok=True)
    result_figure(res, out.with_suffix(".png"), title=f"{Path(image).name}  |  {ref.name}: "
                                                      f"{ref.meta.get('description') or ref.meta.get('stack', '')}")
    rows = summary_rows(res)
    with open(out.with_name(out.name + "_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    maps = dict(tcd=res.tcd, value=res.value, residual=res.residual, alt_value=res.alt_value,
                alt_residual=res.alt_residual)
    if ref.classes:
        maps["classes"] = layer_classes(res)
    np.savez_compressed(out.with_name(out.name + "_maps.npz"), **maps)
    info = dict(image=str(image), reference=str(args.ref), label=ref.label,
                anchors=[a.__dict__ for a in res.anchors], correction_matrix=res.correction.tolist(),
                settings=res.settings)
    out.with_name(out.name + "_run.json").write_text(json.dumps(info, indent=2, default=str))

    print(f"{Path(image).name}: substrate ROI {substrate}, anchors {[(a.name, a.roi) for a in res.anchors[1:]]}")
    print(f"  median residual (assigned) = {np.nanmedian(res.residual[res.assigned]):.2f}, "
          f"unassigned = {100 * (~res.assigned).mean():.1f} %")
    for r in rows:
        if r["fraction"] > 0.005:
            print(f"  {r['item']:>14s}: {100 * r['fraction']:5.1f} %   median residual {r['median_residual']:.2f}")
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


def cmd_materials(args):
    from tcd.materials import library_materials, nk
    from tcd.system import SYSTEMS_DIR, read_system_file

    if args.show:
        wl = np.array([400, 450, 500, 550, 600, 650, 700, 750, 800], float)
        n = nk(args.show, wl)
        print(f"{args.show}:")
        for w, v in zip(wl, n):
            print(f"  {w:5.0f} nm   n = {v.real:.4f}   k = {v.imag:.4f}")
        return
    print("n,k library (data/nk_library.csv), usable by name in system files:")
    print("  " + ", ".join(library_materials()))
    print("\nBundled systems (systems/), usable as --system <name>:")
    for f in sorted(SYSTEMS_DIR.glob("*.jsonc")):
        try:
            desc = read_system_file(f).get("description", "")
        except Exception as e:  # noqa: BLE001 - listing only
            desc = f"(cannot read: {e})"
        print(f"  {f.stem:<16s} {desc}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def system_args(q):
        q.add_argument("--system", required=True, help="system file, or the name of one in systems/ (moo3, ps, ...)")
        q.add_argument("--set", action="append", metavar="NAME=VALUE",
                       help="change a constant of the system file, e.g. --set oxide_nm=285 (repeatable)")
        q.add_argument("--illuminant", help="override optics.illuminant: D65, A, D50 or a CSV lamp spectrum")
        q.add_argument("--na", type=float, help="override optics.na (objective NA; 0 = normal incidence)")

    b = sub.add_parser("build-ref", help="simulate the colour reference of a system file")
    system_args(b)
    b.add_argument("--out", help="reference CSV (default refs/<name>[_<changes>]_<illuminant>.csv)")
    b.add_argument("--spectra", action="store_true", help="also save the reflectance spectra")
    b.set_defaults(func=cmd_build_ref)

    s = sub.add_parser("spectrum", help="reflectance and colour of chosen structures of a system file")
    system_args(s)
    s.add_argument("--at", action="append", metavar="P=V[,Q=W]",
                   help="one structure, e.g. --at thickness_nm=120 (repeatable; default: first, middle, last)")
    s.add_argument("--out", help="figure path (default results/<name>_spectra.png)")
    s.add_argument("--csv", action="store_true", help="also write the spectra as CSV next to the figure")
    s.set_defaults(func=cmd_spectrum)

    m = sub.add_parser("map", help="thickness / layer map of a micrograph")
    m.add_argument("--ref", required=True, help="reference CSV from build-ref")
    m.add_argument("--image", help="micrograph (default: choose in a file dialog)")
    m.add_argument("--out", help="output path stem (default results/<image name>)")
    m.add_argument("--substrate", help="x,y,w,h | auto | colour:r,g,b")
    m.add_argument("--substrate-colour", help="r,g,b (0-1) of the bare substrate; region is found automatically")
    m.add_argument("--anchor", action="append",
                   help="extra calibration region of known structure, e.g. 1L@120,40,30,30, 250nm@colour:r,g,b "
                        "or layers=2,packing=1@x,y,w,h")
    m.add_argument("--interactive", action="store_true", help="click regions that were not given")
    m.add_argument("--gamma", default="1.0",
                   help="camera decoding: 'srgb' or a number g (linear = value**g); default 1.0 = linear camera")
    m.add_argument("--roi-size", type=int, default=40, help="side (px) of auto-placed substrate region")
    m.add_argument("--anchor-size", type=int, default=12, help="side (px) of auto-placed anchor regions")
    m.add_argument("--correction", default="auto", choices=["auto", "diagonal", "matrix"])
    m.add_argument("--max-residual", type=float, default=15.0)
    m.add_argument("--sigma", type=float, default=1.2, help="Gaussian smoothing (px)")
    m.add_argument("--max-side", type=int, default=1600, help="downsample longer side to this (0 = off)")
    m.add_argument("--max-value", "--t-max", "--max-layers", dest="max_value", type=float,
                   help="only consider candidates whose label is <= this (prior knowledge, e.g. from AFM)")
    m.add_argument("--gap", type=float, help="ambiguity gap in label units (default: the reference's)")
    m.set_defaults(func=cmd_map)

    g = sub.add_parser("gamma", help="measure the camera's decoding exponent from an exposure series")
    g.add_argument("--images", nargs="+", required=True, help="same field, only exposure time changed")
    g.add_argument("--exposures", nargs="+", type=float, required=True, help="exposure times, same order")
    g.add_argument("--roi", help="x,y,w,h of an evenly lit region (default: central half of the frame)")
    g.set_defaults(func=cmd_gamma)

    mt = sub.add_parser("materials", help="list the n,k library and the bundled systems")
    mt.add_argument("--show", metavar="NAME", help="print n and k of one library material")
    mt.set_defaults(func=cmd_materials)

    args = p.parse_args(argv)
    if getattr(args, "max_side", None) == 0:
        args.max_side = None
    try:
        args.func(args)
    except (ValueError, KeyError, FileNotFoundError) as e:
        from tcd.expr import ExpressionError
        from tcd.system import SystemError as SysErr
        if isinstance(e, (SysErr, ExpressionError, KeyError, FileNotFoundError)):
            raise SystemExit(f"error: {e}") from None
        raise


if __name__ == "__main__":
    main()
