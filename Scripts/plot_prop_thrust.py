#!/usr/bin/env python3
"""Plot propeller performance vs rotor blade count - thrust, thrust coefficient
and propulsive efficiency - from the swept CFD solutions.

The three quantities are produced by the MATLAB helper make_prop_thrust.m, which
integrates the rotor axial force (momentum + pressure flux across the blade row,
LE -> TE, over the blade span hub -> tip Rc) and the shaft power from the
stagnation-enthalpy rise, then writes <des_base>_thrust_sweep.mat holding

    counts, Tvec, CT, CP, eta, J, Vf, rpm, Ndes, Rc, D, ...

This script reads that file and draws the polished 3-panel figure. It replaces
plot_prop_thrust.m: headless MATLAB (-batch, no display) mangles polylines in
exportgraphics, so the connecting line drops out mid-sweep. matplotlib renders
it reliably and matches the toolchain's Python plotting pattern (plot_prop_chic.py).
Run make_prop_thrust.m first if the *_thrust_sweep.mat is missing.

Usage:
    python plot_prop_thrust.py                          # default design
    python plot_prop_thrust.py 02062026PropPhi06HYBRID
    python plot_prop_thrust.py --dir /path/to/TURBOSTREAM --out /path/perf.pdf
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
    Tvec = np.atleast_1d(d["Tvec"]).astype(float)
    CT = np.atleast_1d(d["CT"]).astype(float)
    eta = np.atleast_1d(d["eta"]).astype(float)
    Vf = float(np.atleast_1d(d["Vf"])[0])
    rpm = float(np.atleast_1d(d["rpm"])[0])
    J = float(np.atleast_1d(d["J"])[0])
    Ndes = int(np.atleast_1d(d["Ndes"])[0])

    panels = [
        (Tvec, r"Thrust,  $T$  [N]", "(a)  thrust"),
        (CT, r"Thrust coefficient,  $C_T$", "(b)  thrust coefficient"),
        (eta, r"Efficiency,  $\eta$", "(c)  efficiency"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
    fig.subplots_adjust(left=0.06, right=0.985, top=0.83, bottom=0.13, wspace=0.28)

    ides = np.where(counts == Ndes)[0]
    for ax, (y, ylab, ttl) in zip(axes, panels):
        ax.plot(counts, y, "-", lw=2.2, color=COL, zorder=2)
        h_cfd, = ax.plot(counts, y, "o", ms=6.5, mfc=COL, mec=COL, zorder=3, label="CFD sweep")
        h_des = None
        if ides.size:
            h_des, = ax.plot(counts[ides], y[ides], "*", ms=15, mfc=CDES, mec="white",
                             mew=0.8, zorder=4, label=f"design (N = {Ndes})")
        ax.grid(True, alpha=0.25, lw=0.7)
        ax.set_axisbelow(True)
        ax.set_xticks(counts)
        ax.set_xlim(counts[0] - 1, counts[-1] + 1)
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(0.0, lo), hi + 0.06 * (hi - lo))
        ax.set_xlabel(r"Rotor blade count,  $N$", fontsize=12.5)
        ax.set_ylabel(ylab, fontsize=13)
        ax.set_title(ttl, fontsize=12.5, fontweight="bold")
        ax.tick_params(labelsize=11)
        if ax is axes[0]:
            handles = [h_cfd] + ([h_des] if h_des is not None else [])
            ax.legend(handles=handles, loc="lower right", fontsize=10, frameon=False)

    fig.suptitle(
        rf"Propeller performance vs blade count"
        rf"        $V_f$ = {Vf:g} m/s,   RPM = {rpm:.0f},   $J$ = {J:.2f}",
        fontsize=14, fontweight="bold")

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropPerf_vs_BladeCount.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
