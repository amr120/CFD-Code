#!/usr/bin/env python3
"""Plot propeller performance vs rotor SOLIDITY - thrust, thrust coefficient and
propulsive efficiency - from the swept CFD solutions.

Same data and quantities as plot_prop_thrust.py, but the abscissa is the rotor
solidity s = C/S instead of the raw blade count N. The blade-count sweep keeps the
blade sections fixed and changes only N, so the pitch S = 2*pi*r_m/N shrinks and
the meanline solidity scales linearly with blade count:

    s(N) = N * C_m / (2*pi*r_m)

The meanline chord C_m is measured from the blade geometry by make_prop_thrust.m,
which writes the solidity vector into <design>_thrust_sweep.mat. (Note: the mission
JSON's 'sigma' is the Chapter-3 disk-to-jet velocity ratio, NOT the solidity.)
Blade-count tick labels are kept on a secondary top axis so both metrics are read.

Usage:
    python plot_prop_solidity.py                          # default design
    python plot_prop_solidity.py 02062026PropPhi06HYBRID
    python plot_prop_solidity.py --dir /path/to/TURBOSTREAM --out /path/perf.pdf
"""

import argparse
import os

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# 14 pt base font everywhere and a thicker black axes frame (cf. plot_prop_chic.py)
plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_OUT_DIR = "/Data/Engine_Selector"
DEFAULT_DESIGN = "02062026PropPhi06HYBRID"

COL = "#2059a8"        # sweep line / markers
CDES = "#d9541a"       # design-point highlight
CNEW = "#555555"       # Newman s = 0.5 boundary

NEWMAN_S = 0.5         # Newman (2021) propeller/fan solidity convention
GOLDEN = (1.0 + 5.0 ** 0.5) / 2.0        # each plot panel is a golden rectangle


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design", nargs="?", default=DEFAULT_DESIGN,
                    help="design base name (default: %(default)s)")
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM dir with the .mat")
    ap.add_argument("--out", default=None, help="output path (png; a matching pdf is also written)")
    args = ap.parse_args()

    mat = os.path.join(args.dir, args.design + "_thrust_sweep.mat")
    d = sio.loadmat(mat, squeeze_me=True)
    counts = np.atleast_1d(d["counts"]).astype(float)
    sol = np.atleast_1d(d["solidity"]).astype(float)
    Tvec = np.atleast_1d(d["Tvec"]).astype(float)
    CT = np.atleast_1d(d["CT"]).astype(float)
    eta = np.atleast_1d(d["eta"]).astype(float)
    Vf = float(np.atleast_1d(d["Vf"])[0])
    rpm = float(np.atleast_1d(d["rpm"])[0])
    J = float(np.atleast_1d(d["J"])[0])
    Ndes = int(np.atleast_1d(d["Ndes"])[0])

    panels = [
        (Tvec, r"Thrust,  $T$  [N]"),
        (CT, r"Thrust coefficient,  $C_T$"),
        (eta, r"Efficiency,  $\eta$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 13.0 / GOLDEN))   # overall box = golden
    fig.subplots_adjust(left=0.06, right=0.985, top=0.88, bottom=0.09, wspace=0.28)

    ides = np.where(counts == Ndes)[0]
    for ax, (y, ylab) in zip(axes, panels):
        # Newman s = 0.5 propeller/fan boundary
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
        ax.set_xlim(sol[0] - dx, sol[-1] + dx)
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(0.0, lo), hi + 0.06 * (hi - lo))
        ax.set_xlabel(r"Rotor solidity,  $s = C/S$", fontsize=12.5)
        ax.set_ylabel(ylab, fontsize=13)
        ax.tick_params(labelsize=11)

        # label Newman's boundary at the top of the axes
        ytop = ax.get_ylim()[1]
        ax.annotate(r"Newman $s=0.5$", xy=(NEWMAN_S, ytop), xytext=(4, -3),
                    textcoords="offset points", rotation=90, va="top", ha="left",
                    fontsize=8.5, color=CNEW)

        # secondary top axis: blade-count labels at the same solidity positions
        axt = ax.secondary_xaxis("top")
        axt.set_xticks(sol)
        axt.set_xticklabels([f"{int(c)}" for c in counts], fontsize=9.5)
        axt.set_xlabel(r"blade count,  $N$", fontsize=10, labelpad=3)
        axt.tick_params(length=3)

        if ax is axes[0]:
            handles = [h_cfd] + ([h_des] if h_des is not None else [])
            ax.legend(handles=handles, loc="lower right", fontsize=9.5, frameon=False)

    fig.suptitle(
        "Propeller performance vs rotor solidity",
        fontsize=14, fontweight="bold", y=0.98)

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropPerf_vs_Solidity.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
