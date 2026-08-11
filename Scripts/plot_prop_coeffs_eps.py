#!/usr/bin/env python3
"""Plot the measured performance coefficients against the flow-contraction
coefficient epsilon - the three-dimensional-rotor analogue of the Chapter-3
epsilon sweep (Fig.~\\ref{fig:epsilon-sweep}), which plots the ideal coefficients
against epsilon.

For each blade count of the sweep the measured coefficients
    C_T   = Th / (mdot * V_j)        thrust coefficient   (make_prop_eps: CT_new)
    C_P   = P  / (mdot * V_j^2)      power coefficient    (make_prop_eps: CP_new)
    eta   = Th * V_fs / P            propulsive efficiency (make_prop_thrust: eta)
are plotted against the measured epsilon = V_fs / V_j. In the unified framework
these are single-valued functions of epsilon, so the swept points collapse onto
single curves; here they are the measured (viscous, three-dimensional) values,
not the ideal actuator-disk relations. All are read from
<design>_thrust_sweep.mat (written by make_prop_thrust.m / make_prop_eps.m).

Usage:
    python plot_prop_coeffs_eps.py
    python plot_prop_coeffs_eps.py 02062026PropPhi06HYBRID
    python plot_prop_coeffs_eps.py --dir /path/to/TURBOSTREAM --out /path/fig.pdf
"""

import argparse
import os

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_OUT_DIR = "/Data/Engine_Selector"
DEFAULT_DESIGN = "02062026PropPhi06HYBRID"

GOLDEN = (1.0 + 5.0 ** 0.5) / 2.0        # overall figure box = golden ratio

# One series per curve: mat key, label, colour, marker, ideal Ch-3 curve f(eps)
SERIES = [
    ("CT_new", r"$C_T$",                        "#1f77b4", "o", lambda e: 1.0 - e),
    ("CP_new", r"$C_P$",                        "#2ca02c", "s", lambda e: 0.5 * (1.0 - e ** 2)),
    ("eta",    r"$\eta_{\mathrm{propulsive}}$", "#d62728", "D", lambda e: 2.0 / (1.0 + 1.0 / e)),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design", nargs="?", default=DEFAULT_DESIGN)
    ap.add_argument("--dir", default=DEFAULT_TS_DIR)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = sio.loadmat(os.path.join(args.dir, args.design + "_thrust_sweep.mat"),
                    squeeze_me=True)
    eps = np.atleast_1d(d["eps"]).astype(float)
    order = np.argsort(eps)                       # plot lines in increasing epsilon

    fig, ax = plt.subplots(figsize=(8.2, 8.2 / GOLDEN))    # overall box = golden ratio
    fig.subplots_adjust(left=0.10, right=0.965, top=0.965, bottom=0.13)

    ee = np.linspace(eps.min(), eps.max(), 200)   # smooth eps for the ideal curves
    for key, lab, col, mk, ideal in SERIES:
        y = np.atleast_1d(d[key]).astype(float)
        ax.plot(ee, ideal(ee), "--", lw=1.6, color=col, alpha=0.9, zorder=1.5)  # ideal (Ch.3)
        ax.plot(eps[order], y[order], "-", lw=2.0, color=col, zorder=2)         # measured
        ax.plot(eps, y, mk, ms=5.5, mfc=col, mec=col, zorder=3, label=lab)

    ax.grid(True, alpha=0.25, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_ylim(0.0, 0.75)
    ax.set_xlabel(r"Flow contraction,  $\varepsilon = V_{\mathrm{fs}}/V_{\mathrm{j}}$",
                  fontsize=14)
    ax.set_ylabel(r"[ - ]", fontsize=14)
    ax.tick_params(labelsize=14)

    # legend: the three coefficients (colour) plus the measured/ideal line styles
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([0], [0], color="0.35", ls="-", lw=2.0),
                Line2D([0], [0], color="0.35", ls="--", lw=1.6)]
    labels += ["measured", "ideal (Ch. 3)"]
    ax.legend(handles, labels, loc="lower right", fontsize=13, handlelength=1.8,
              ncol=2, columnspacing=1.3, frameon=True, facecolor="white",
              framealpha=0.85, edgecolor="0.8")

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropCoeffs_vs_Eps.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
