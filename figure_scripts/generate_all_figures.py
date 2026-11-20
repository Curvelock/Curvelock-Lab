"""
generate_all_figures.py
=======================
Driver for the Paper 1 public figure package.

Behaviour:
  * creates the output directory
  * runs each figure generator in order
  * catches per-figure exceptions and continues (so one missing data product
    does not block the others)
  * prints a summary table of generated / failed figures
  * exits nonzero if any MANDATORY figure failed
  * NEVER generates placeholder or dummy figures for missing data - a missing
    figure is reported missing

Figure 1 is mandatory-and-dataless (it always generates). Figures 2-7 are
mandatory but data-dependent; if their data products are absent, they fail
loudly and the driver reports them missing.

Usage:
    python generate_all_figures.py --data data --out figures_output
    python generate_all_figures.py --data <dir> --out <dir> --synthetic
        (the --synthetic flag only routes the watermark; it does not
         fabricate data and is intended for the development self-test)
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import traceback

import matplotlib
matplotlib.use("Agg")
import common_style as cs
import figure1_schematic
import figure2_spectrum
import figure3_modes
import figure4_heatmap
import figure5_capture
import figure6_budget
import figure7_scan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=cs.DATA_ROOT,
                    help="path to the data/ contract root")
    ap.add_argument("--out", type=Path, default=cs.OUT_ROOT,
                    help="output directory (created if absent)")
    ap.add_argument("--synthetic", action="store_true",
                    help="watermark figures as synthetic (dev self-test only)")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    cs.write_manifest(args.out, [])
    syn = args.synthetic

    # (name, callable, mandatory)
    jobs = [
        ("figure1_schematic", lambda: figure1_schematic.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure2_spectrum", lambda: figure2_spectrum.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure3_modes", lambda: figure3_modes.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure4_heatmap", lambda: figure4_heatmap.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure5_capture", lambda: figure5_capture.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure6_budget", lambda: figure6_budget.make_figure(
            args.data, args.out, synthetic=syn), True),
        ("figure7_scan", lambda: figure7_scan.make_figure(
            args.data, args.out, synthetic=syn), True),
    ]

    results = []
    for name, fn, mandatory in jobs:
        try:
            files = fn()
            results.append((name, "OK", mandatory, files, None))
        except Exception as exc:  # keep going; report at the end
            results.append((name, "FAILED", mandatory, [],
                            f"{type(exc).__name__}: {exc}"))
            traceback.print_exc()

    # ---- summary table -----------------------------------------------
    print("\n" + "=" * 68)
    print(f"{'figure':<22}{'status':<10}{'mandatory':<11}outputs")
    print("-" * 68)
    any_mandatory_failed = False
    for name, status, mandatory, files, err in results:
        out = ", ".join(files) if files else (err or "")
        print(f"{name:<22}{status:<10}{str(mandatory):<11}{out}")
        if status == "FAILED" and mandatory:
            any_mandatory_failed = True
    print("=" * 68)
    if syn:
        print("NOTE: figures were watermarked SYNTHETIC (dev self-test). "
              "These are NOT scientific results.")
    print(f"manifest: {args.out / 'figures_manifest.json'}")

    if any_mandatory_failed:
        print("\nOne or more MANDATORY figures failed. "
              "Ensure the public simulation data products exist under "
              f"'{args.data}'. No figures were faked.")
        sys.exit(1)
    print("\nAll mandatory figures generated.")


if __name__ == "__main__":
    main()
