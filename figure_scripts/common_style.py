"""Shared plotting, output paths and provenance for Paper 1 figures."""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = ROOT / "data"
OUT_ROOT = ROOT / "figures_output"
PACKAGE_VERSION = "1.2-integrated"
PALETTE = {"total": "#000000", "analytic": "#555555", "well_fill": "#dcdcdc"}
MODE_COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
MODE_STYLES = ("-", "--", "-.", ":")


def isolated_style(fn):
    """Do not let an individual generator's rcParams leak into another."""
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        before = set(plt.get_fignums())
        with mpl.rc_context({"font.family": "serif", "mathtext.fontset": "stix",
                             "font.size": 9, "axes.labelsize": 10,
                             "xtick.labelsize": 9, "ytick.labelsize": 9,
                             "legend.fontsize": 8, "xtick.direction": "in",
                             "ytick.direction": "in", "xtick.top": True,
                             "ytick.right": True, "axes.linewidth": 0.9,
                             "lines.linewidth": 1.3}):
            try:
                return fn(*args, **kwargs)
            finally:
                for num in set(plt.get_fignums()) - before:
                    plt.close(num)
    return wrapped


def figsize(columns, height):
    return (3.4 if columns == 1 else 7.0, height)


def mode_color(n):
    return MODE_COLORS[(n - 1) % len(MODE_COLORS)]


def mode_linestyle(n):
    return MODE_STYLES[(n - 1) % len(MODE_STYLES)]


def shade_well(ax, a, ymin=0, ymax=1):
    ax.axvspan(0, a, ymin=ymin, ymax=ymax, color=PALETTE["well_fill"], zorder=0)
    ax.axvline(a, color=PALETTE["analytic"], ls="--", lw=0.8, zorder=1)


def gate_line(ax, y, label):
    color = "#D55E00"
    ax.axhline(y, color=color, ls="--", lw=0.9)
    ax.text(1.005, y, label, color=color, va="center", ha="left",
            transform=ax.get_yaxis_transform(), clip_on=False)


def write_manifest(out_dir, records):
    path = Path(out_dir) / "figures_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")


def save_figure(fig, name, out_dir, source_data, run_ids=(), synthetic=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if synthetic:
        fig.text(0.5, 0.5, "SYNTHETIC — NOT SCIENTIFIC DATA", ha="center",
                 va="center", color="crimson", alpha=0.45, rotation=25)
    files = [name + ".pdf", name + ".png"]
    for filename in files:
        fig.savefig(out / filename, dpi=600, bbox_inches="tight", pad_inches=0.05)
    sources = []
    for path in map(Path, source_data):
        try:
            label = path.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            label = path.resolve().as_posix()
        sources.append({"path": label, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    manifest = out / "figures_manifest.json"
    records = json.loads(manifest.read_text(encoding="utf-8")) if manifest.exists() else []
    records = [r for r in records if r["figure"] != name]
    records.append({"figure": name, "files": files,
                    "source_data": [source["path"] for source in sources],
                    "source_sha256": {source["path"]: source["sha256"] for source in sources},
                    "run_ids": [r for r in run_ids if r is not None],
                    "figure_package_version": PACKAGE_VERSION,
                    "synthetic": bool(synthetic)})
    write_manifest(out, sorted(records, key=lambda r: r["figure"]))
    return [str(out / filename) for filename in files]


def standalone(make_figure):
    parser = argparse.ArgumentParser(description=make_figure.__module__)
    parser.add_argument("--data", type=Path, default=DATA_ROOT)
    parser.add_argument("--out", type=Path, default=OUT_ROOT)
    parser.add_argument("--synthetic", action="store_true")
    args = parser.parse_args()
    print("\n".join(make_figure(args.data, args.out, synthetic=args.synthetic)))
