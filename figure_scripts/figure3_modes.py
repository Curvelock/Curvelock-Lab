"""
figure3_modes.py
================
Figure 3 - Localized standing-wave eigenmodes of the prescribed well.

Panel A : normalized mode profiles u_n(r), well region shaded.
Panel B : exterior decay log10|u_n(r)| for r > a, with dashed reference
          slopes -kappa_n/ln(10) demonstrating exponential localization,
          kappa_n = sqrt(m0^2 - omega_n^2).

Only ordinary bound modes (0 < omega_n^2 < m0^2) are shown.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import common_style as cs
import common_io as cio


@cs.isolated_style
def make_figure(data_dir, out_dir, m0=1.0, floor_exp=-8.0, synthetic=False):
    base = f"{data_dir}/baseline"
    prof_path = f"{base}/mode_profiles.csv"
    meta_path = f"{base}/run_metadata.json"

    df, u_cols = cio.load_mode_profiles(prof_path)
    meta = cio.load_run_metadata(meta_path)
    modes = {int(m["n"]): m for m in meta["modes"]}

    r = df["r"].to_numpy()
    a = float(meta["params"]["a"]) if "params" in meta else float(meta["a"])

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=cs.figsize(1, height=4.2))

    # ---- Panel A: profiles --------------------------------------------
    cs.shade_well(axA, a)
    for col in u_cols:
        n = int(col.split("_")[-1])
        axA.plot(r, df[col], ls=cs.mode_linestyle(n), color=cs.mode_color(n),
                 label=f"mode {n}")
    axA.axhline(0.0, color="#999999", lw=0.5, zorder=0)
    axA.set_ylabel(r"$u_n(r)$  (normalized)")
    axA.set_xlim(0, min(4.0 * a, r.max()))
    axA.legend(loc="upper right", frameon=False)

    # ---- Panel B: exterior decay --------------------------------------
    cs.shade_well(axB, a, ymin=0, ymax=1)
    ext = r >= a
    for col in u_cols:
        n = int(col.split("_")[-1])
        u = np.abs(df[col].to_numpy())
        logu = np.log10(np.clip(u, 10.0 ** floor_exp, None))
        axB.plot(r[ext], logu[ext], ls=cs.mode_linestyle(n),
                 color=cs.mode_color(n), label=f"mode {n}")
        # analytic decay slope reference, if omega2 available
        m = modes.get(n)
        if m is not None and "omega2" in m:
            kappa = np.sqrt(max(m0 ** 2 - float(m["omega2"]), 0.0))
            # anchor the reference line at r = a to the mode's value there
            i_a = np.argmin(np.abs(r - a))
            anchor = np.log10(max(u[i_a], 10.0 ** floor_exp))
            ref = anchor - (kappa / np.log(10.0)) * (r[ext] - a)
            axB.plot(r[ext], ref, ls=(0, (2, 2)), lw=0.8,
                     color=cs.PALETTE["analytic"], zorder=1)
    axB.text(0.98, 0.95, r"dashed: $e^{-\kappa_n r}$,  "
             r"$\kappa_n=\sqrt{m_0^2-\omega_n^2}$",
             transform=axB.transAxes, ha="right", va="top", fontsize=7,
             color=cs.PALETTE["analytic"])
    axB.set_xlabel(r"$r$  [model units]")
    axB.set_ylabel(r"$\log_{10}|u_n(r)|$")
    axB.set_xlim(a, r.max())
    axB.set_ylim(floor_exp, 0.4)

    fig.align_ylabels([axA, axB])
    return cs.save_figure(
        fig, "figure3_modes", out_dir,
        source_data=[prof_path, meta_path],
        run_ids=[meta.get("run_id")], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)
