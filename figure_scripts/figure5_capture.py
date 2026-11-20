"""
figure5_capture.py
==================
Figure 5 - Nonlinear loading of the bound-mode subspace.

f_cap(t) (total) and per-mode contributions E_n(t)/E_in for the baseline
nonlinear run, with the SAME diagnostic for the linear null run (lambda_chi=0)
overlaid as the numerical floor. Gate lines eps_lin and f_min are read from
metadata. The final analysis window is shaded.

The caption must cross-reference Figure 6: capture is adjudicated jointly with
energy-budget closure, never by f_cap growth alone.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import common_style as cs
import common_io as cio


@cs.isolated_style
def make_figure(data_dir, out_dir, y_floor=1e-6, synthetic=False):
    base = f"{data_dir}/baseline"
    nullp = f"{data_dir}/null_linear"
    ts_path = f"{base}/time_series.csv"
    meta_path = f"{base}/run_metadata.json"
    null_ts_path = f"{nullp}/time_series.csv"

    ts = cio.load_time_series(ts_path)
    meta = cio.load_run_metadata(meta_path)
    gates = cio.get_gates(meta)

    # linear null run is mandatory for this figure
    null_ts = cio.load_time_series(null_ts_path)

    fig, ax = plt.subplots(figsize=cs.figsize(1, height=2.9))

    t = ts["t"].to_numpy()

    # total f_cap
    ax.semilogy(t, np.clip(ts["f_cap"], y_floor, None),
                color=cs.PALETTE["total"], lw=1.3, label=r"$f_{\mathrm{cap}}(t)$")

    # per-mode contributions
    E_in_cols = cio.resolved_mode_indices(ts)
    for n in E_in_cols:
        # E_mode_n columns are absolute energies; f_cap is the summed fraction.
        # Per-mode fraction = E_mode_n / E_in, and E_in = sum_n E_mode_n / f_cap.
        frac = ts[f"E_mode_{n}"].to_numpy()
        denom = np.where(ts["f_cap"].to_numpy() > 0,
                         np.nansum([ts[f"E_mode_{m}"].to_numpy()
                                    for m in E_in_cols], axis=0)
                         / np.clip(ts["f_cap"].to_numpy(), 1e-300, None),
                         np.nan)
        mode_frac = np.where(denom > 0, frac / denom, np.nan)
        ax.semilogy(t, np.clip(mode_frac, y_floor, None),
                    ls=cs.mode_linestyle(n), color=cs.mode_color(n),
                    lw=1.0, label=f"mode {n}")

    # linear-null floor
    ax.semilogy(null_ts["t"], np.clip(null_ts["f_cap"], y_floor, None),
                color="#888888", lw=0.9, ls=(0, (1, 1)),
                label=r"linear null ($\lambda_\chi=0$)")

    # gates
    cs.gate_line(ax, gates["eps_lin"], r"$\varepsilon_{\mathrm{lin}}$")
    cs.gate_line(ax, gates["f_min"], r"$f_{\min}$")

    # final analysis window
    win = meta["outcomes"].get("analysis_window")
    if win and len(win) == 2:
        ax.axvspan(win[0], win[1], color=cs.PALETTE["well_fill"], zorder=0)

    ax.set_xlabel(r"$t$  [model units]")
    ax.set_ylabel(r"fraction of $E_{\mathrm{in}}$")
    ax.set_ylim(y_floor, 1.2)
    ax.set_xlim(t.min(), t.max())
    ax.legend(loc="lower right", frameon=False, ncol=1, fontsize=7)

    return cs.save_figure(
        fig, "figure5_capture", out_dir,
        source_data=[ts_path, null_ts_path, meta_path],
        run_ids=[meta.get("run_id")], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)
