#!/usr/bin/env python3
"""Plot the Chapter-3 thrust coefficient and the flow-contraction coefficient
against rotor SOLIDITY, as a two-panel figure.

Left panel: the reformulated thrust coefficient of Chapter 3,
    C_T = Thrust / (mdot * V_j),
computed FRESHLY from the CFD - the integrated rotor thrust over the rotor mass
flow and the measured developed-jet velocity (make_prop_eps.m writes CT_new). It
is NOT evaluated as 1 - eps; it equals 1 - eps only to the extent the momentum
relation Th = mdot*(V_j - V_fs) holds, so the small offset from 1 - eps is the
residual of that relation.

Right panel: the flow-contraction coefficient eps = V_fs / V_j (the same quantity
as plot_prop_eps.py), the Chapter-3 loading ratio.

Both are drawn against the meanline solidity s = C/S, with the blade count on a
secondary axis and Newman's s = 0.5 convention marked, in the golden-ratio,
title-clean style of the other Section 6.2 figures.

Usage:
    python plot_prop_ct_eps.py
    python plot_prop_ct_eps.py 02062026PropPhi06HYBRID
    python plot_prop_ct_eps.py --dir /path/to/TURBOSTREAM --out /path/fig.pdf
"""

import argparse
import os

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_OUT_DIR = "/Data/Engine_Selector"
DEFAULT_DESIGN = "02062026PropPhi06HYBRID"

COL = "#2059a8"
CDES = "#d9541a"
CNEW = "#555555"
NEWMAN_S = 0.5
GOLDEN = (1.0 + 5.0 ** 0.5) / 2.0        # overall figure box = golden ratio


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design", nargs="?", default=DEFAULT_DESIGN)
    ap.add_argument("--dir", default=DEFAULT_TS_DIR)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = sio.loadmat(os.path.join(args.dir, args.design + "_thrust_sweep.mat"),
                    squeeze_me=True)
    counts = np.atleast_1d(d["counts"]).astype(float)
    sol = np.atleast_1d(d["solidity"]).astype(float)
    CT = np.atleast_1d(d["CT_new"]).astype(float)
    eps = np.atleast_1d(d["eps"]).astype(float)
    Ndes = int(np.atleast_1d(d["Ndes"])[0])

    # keep counts with finite values (in case any solution was missing)
    ok = np.isfinite(CT) & np.isfinite(eps)
    counts, sol, CT, eps = counts[ok], sol[ok], CT[ok], eps[ok]

    panels = [
        (CT, r"Thrust coefficient,  $C_T = \mathrm{Th}/(\dot{m}\,V_{\mathrm{j}})$"),
        (eps, r"Flow contraction,  $\varepsilon = V_{\mathrm{fs}}/V_{\mathrm{j}}$"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 10.6 / GOLDEN))   # overall box = golden
    fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.10, wspace=0.24)

    ides = np.where(counts == Ndes)[0]
    for ax, (y, ylab) in zip(axes, panels):
        ax.axvline(NEWMAN_S, color=CNEW, ls="--", lw=1.4, zorder=1)

        ax.plot(sol, y, "-", lw=2.2, color=COL, zorder=2)
        h_cfd, = ax.plot(sol, y, "o", ms=6.5, mfc=COL, mec=COL, zorder=3, label="CFD sweep")
        h_des = None
        if ides.size:
            h_des, = ax.plot(sol[ides], y[ides], "*", ms=15, mfc=CDES, mec="white",
                             mew=0.8, zorder=4,
                             label=rf"design ($s$ = {sol[ides][0]:.2f}, N = {Ndes})")
        ax.grid(True, alpha=0.25, lw=0.7)
        ax.set_axisbelow(True)
        dx = 0.5 * (sol[1] - sol[0])
        ax.set_xlim(min(sol[0] - dx, NEWMAN_S - 0.1), sol[-1] + dx)
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(0.0, lo), hi + 0.08 * (hi - lo))
        ax.set_xlabel(r"Rotor solidity,  $s = C/S$", fontsize=12.5)
        ax.set_ylabel(ylab, fontsize=12.5)
        ax.tick_params(labelsize=11)

        ytop = ax.get_ylim()[1]
        ax.annotate(r"Newman $s=0.5$", xy=(NEWMAN_S, ytop), xytext=(4, -3),
                    textcoords="offset points", rotation=90, va="top", ha="left",
                    fontsize=9, color=CNEW)

        axt = ax.secondary_xaxis("top")
        axt.set_xticks(sol)
        axt.set_xticklabels([f"{int(c)}" for c in counts], fontsize=9.5)
        axt.set_xlabel(r"blade count,  $N$", fontsize=10, labelpad=3)
        axt.tick_params(length=3)

        if ax is axes[0]:
            handles = [h_cfd] + ([h_des] if h_des is not None else [])
            ax.legend(handles=handles, loc="lower right", fontsize=9.5, frameon=False)

    fig.suptitle(
        "3D rotor thrust coefficient and flow contraction vs solidity",
        fontsize=13.5, fontweight="bold", y=0.98)

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropCT_Eps_vs_Solidity.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
