#!/usr/bin/env python3
"""Plot 3D-rotor flow turning and stagnation-pressure loss against rotor SOLIDITY,
to sit directly alongside the 2D MISES cascade results of thesis Chapter 6.

Chapter 6 established, at the 2D blade-section level, that raising solidity
increases the flow turning and the profile loss smoothly and continuously, with
no feature at Newman's s = 0.5 propeller/fan boundary. This figure is the 3D-rotor
analogue for the HYBRID propeller blade-count sweep (thesis Section 6.2): the
turning and loss are integrated across the real rotor from the swept CFD.

Both quantities are computed in make_prop_thrust.m over the blade span (hub -> tip
Rc), mass-averaged in the rotor-relative frame across the rotor LE and TE planes:

    turning  Delta-beta = | alpha_rel,1 - alpha_rel,2 |            [deg]
    loss     Yp = (Po_rel,1 - Po_rel,2) / (Po_rel,1 - P_1)        (cf. Ch.5 Eq. Yp)

and written, with the solidity vector, into <design>_thrust_sweep.mat. Newman's
s = 0.5 boundary is drawn for direct comparison with the Chapter 6 figures.

Usage:
    python plot_prop_turnloss.py                          # default design
    python plot_prop_turnloss.py 02062026PropPhi06HYBRID
    python plot_prop_turnloss.py --dir /path/to/TURBOSTREAM --out /path/fig.pdf
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
    turn = np.atleast_1d(d["turn"]).astype(float)
    Yp = np.atleast_1d(d["Yp"]).astype(float)
    Vf = float(np.atleast_1d(d["Vf"])[0])
    rpm = float(np.atleast_1d(d["rpm"])[0])
    J = float(np.atleast_1d(d["J"])[0])
    Ndes = int(np.atleast_1d(d["Ndes"])[0])

    panels = [
        (turn, r"Rotor flow turning,  $\Delta\beta$  [deg]"),
        (Yp, r"Loss coefficient,  $Y_p$"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 10.6 / GOLDEN))   # overall box = golden
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.10, wspace=0.26)

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
        ax.set_xlim(min(sol[0] - dx, NEWMAN_S - 0.1), sol[-1] + dx)
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(0.0, lo), hi + 0.08 * (hi - lo))
        ax.set_xlabel(r"Rotor solidity,  $s = C/S$", fontsize=12.5)
        ax.set_ylabel(ylab, fontsize=13)
        ax.tick_params(labelsize=11)

        # label Newman's boundary at the top of the axes
        ytop = ax.get_ylim()[1]
        ax.annotate(r"Newman $s=0.5$", xy=(NEWMAN_S, ytop), xytext=(4, -3),
                    textcoords="offset points", rotation=90, va="top", ha="left",
                    fontsize=9, color=CNEW)

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
        "3D rotor turning and loss vs solidity",
        fontsize=14, fontweight="bold", y=0.98)

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropTurnLoss_vs_Solidity.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
