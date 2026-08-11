#!/usr/bin/env python3
"""Plot the flow-contraction coefficient eps against rotor SOLIDITY for the
blade-count sweep.

eps = V_fs / V_j is the Chapter-3 propulsor loading ratio: freestream velocity
over developed-jet velocity, bounded 0-1, with lower eps meaning a more heavily
loaded propulsor. It is measured from each CFD solution by make_prop_eps.m (the
domain-inlet freestream over the peak capture-streamtube jet velocity behind the
rotor) and appended to <design>_thrust_sweep.mat.

As blades are added the rotor is loaded more heavily, the jet is accelerated and
eps falls - the direct link between solidity and loading that underpins the
propeller/fan spectrum of Chapter 3. This script draws eps vs the meanline
solidity s = C/S, with the blade count on a secondary axis and Newman's s = 0.5
convention marked, in the same style as plot_prop_solidity.py.

Usage:
    python plot_prop_eps.py
    python plot_prop_eps.py 02062026PropPhi06HYBRID
    python plot_prop_eps.py --dir /path/to/TURBOSTREAM --out /path/eps.pdf
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
GOLDEN = (1.0 + 5.0 ** 0.5) / 2.0        # each plot panel is a golden rectangle


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
    eps = np.atleast_1d(d["eps"]).astype(float)
    Ndes = int(np.atleast_1d(d["Ndes"])[0])
    Vf = float(np.atleast_1d(d["Vf"])[0])
    rpm = float(np.atleast_1d(d["rpm"])[0])
    J = float(np.atleast_1d(d["J"])[0])

    # keep only counts with a finite eps (in case any solution was missing)
    ok = np.isfinite(eps)
    counts, sol, eps = counts[ok], sol[ok], eps[ok]

    fig, ax = plt.subplots(figsize=(7.4, 7.4 / GOLDEN))    # overall box = golden ratio
    fig.subplots_adjust(left=0.12, right=0.965, top=0.84, bottom=0.14)

    ax.axvline(NEWMAN_S, color=CNEW, ls="--", lw=1.4, zorder=1)
    ax.plot(sol, eps, "-", lw=2.2, color=COL, zorder=2)
    h_cfd, = ax.plot(sol, eps, "o", ms=6.5, mfc=COL, mec=COL, zorder=3, label="CFD sweep")

    ides = np.where(counts == Ndes)[0]
    handles = [h_cfd]
    if ides.size:
        h_des, = ax.plot(sol[ides], eps[ides], "*", ms=16, mfc=CDES, mec="white",
                         mew=0.8, zorder=4,
                         label=rf"design ($s$ = {sol[ides][0]:.2f}, N = {Ndes})")
        handles.append(h_des)

    ax.grid(True, alpha=0.25, lw=0.7)
    ax.set_axisbelow(True)
    dx = 0.5 * (sol[1] - sol[0])
    ax.set_xlim(min(sol[0] - dx, NEWMAN_S - 0.1), sol[-1] + dx)
    lo, hi = ax.get_ylim()
    ax.set_ylim(max(0.0, lo - 0.03), min(1.0, hi + 0.05))
    ax.set_xlabel(r"Rotor solidity,  $s = C/S$", fontsize=13)
    ax.set_ylabel(r"Flow contraction,  $\varepsilon = V_{\mathrm{fs}}/V_{\mathrm{j}}$",
                  fontsize=13)
    ax.tick_params(labelsize=11)

    ytop = ax.get_ylim()[1]
    ax.annotate(r"Newman $s=0.5$", xy=(NEWMAN_S, ytop), xytext=(4, -4),
                textcoords="offset points", rotation=90, va="top", ha="left",
                fontsize=9, color=CNEW)

    axt = ax.secondary_xaxis("top")
    axt.set_xticks(sol)
    axt.set_xticklabels([f"{int(c)}" for c in counts], fontsize=9)
    axt.set_xlabel(r"blade count,  $N$", fontsize=10, labelpad=3)
    axt.tick_params(length=3)

    ax.legend(handles=handles, loc="upper right", fontsize=10, frameon=False)
    fig.suptitle("Propeller flow contraction vs solidity",
                 fontsize=13, fontweight="bold", y=0.985)

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropEps_vs_Solidity.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
