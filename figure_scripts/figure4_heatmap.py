"""
figure4_heatmap.py
==================
Figure 4 - Space-time evolution log10|chi(t, r)| of the orthogonalized packet.

Shows the ingoing packet, interaction with the well, outgoing radiation
absorbed in the sponge region, and persistent localized oscillation near the
well. The late-time localized structure is labeled only as "persistent
localized oscillation" - never as a particle/object/soliton/mass.

Colormap: perceptually uniform sequential (magma). Never rainbow/jet.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import common_style as cs
import common_io as cio


@cs.isolated_style
def make_figure(data_dir, out_dir, floor_exp=-8.0, synthetic=False):
    base = f"{data_dir}/baseline"
    npz_path = f"{base}/field_evolution.npz"
    meta_path = f"{base}/run_metadata.json"

    r, t, chi = cio.load_field_evolution(npz_path)
    meta = cio.load_run_metadata(meta_path)

    params = meta.get("params", meta)
    a = float(params["a"])
    grid = meta.get("grid", params)
    r_max = float(grid.get("r_max", r.max()))
    sponge = meta.get("sponge", {})
    start_frac = float(sponge.get("start_frac", 0.75))
    r_sponge = start_frac * r_max

    logmag = np.log10(np.clip(np.abs(chi), 10.0 ** floor_exp, None))
    vmax = float(logmag.max())

    fig, ax = plt.subplots(figsize=cs.figsize(2, height=3.6))

    # nearest-neighbour shading: do not manufacture features by interpolation
    mesh = ax.pcolormesh(r, t, logmag, shading="nearest",
                         cmap="magma", vmin=floor_exp, vmax=vmax,
                         rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}|\chi(t,r)|$")

    # well edge + sponge band
    ax.axvline(a, color="white", ls=(0, (4, 3)), lw=0.9)
    ax.text(a, t.max() * 1.005, "well edge", color="#222222", fontsize=7,
            ha="center", va="bottom")
    ax.axvspan(r_sponge, r_max, facecolor="none", edgecolor="white",
               hatch="//", lw=0.0, alpha=0.6)
    ax.text((r_sponge + r_max) / 2, t.max() * 0.02, "absorbing layer",
            color="white", fontsize=7, ha="center", va="bottom")

    # three neutral, approved-vocabulary annotations
    ax.set_xlabel(r"$r$  [model units]")
    ax.set_ylabel(r"$t$  [model units]")
    ax.set_xlim(0, r_max)
    ax.set_ylim(t.min(), t.max())

    return cs.save_figure(
        fig, "figure4_heatmap", out_dir,
        source_data=[npz_path, meta_path],
        run_ids=[meta.get("run_id")], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)
