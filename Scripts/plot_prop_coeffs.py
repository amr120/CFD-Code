#!/usr/bin/env python3
"""Plot the Chapter-3 unified performance coefficients against rotor SOLIDITY,
all on a single axes - the three-dimensional-rotor analogue of the epsilon sweep
of Chapter 3 (Fig.~\\ref{fig:epsilon-sweep}).

Four measured quantities are shown vs the meanline solidity s = C/S:
    C_T   = Th / (mdot * V_j)          thrust coefficient   (make_prop_eps: CT_new)
    C_P   = Power / (mdot * V_j^2)     power coefficient    (make_prop_eps: CP_new)
    eta   = Th * V_fs / Power          propulsive efficiency (make_prop_thrust: eta)
    eps   = V_fs / V_j                 flow contraction     (make_prop_eps: eps)

C_T and C_P are computed FRESHLY from the integrated rotor thrust / shaft power,
the rotor mass flow and the measured developed-jet velocity - not from the ideal
actuator-disk identities C_T = 1 - eps, C_P = 0.5(1 - eps^2). All four are read
from <design>_thrust_sweep.mat.

Newman's s = 0.5 convention and the design point are marked. Single golden-ratio
axes, in the title-clean style of the other Section 6.2 figures.

Usage:
    python plot_prop_coeffs.py
    python plot_prop_coeffs.py 02062026PropPhi06HYBRID
    python plot_prop_coeffs.py --dir /path/to/TURBOSTREAM --out /path/fig.pdf
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

CNEW = "#555555"                         # Newman s = 0.5 boundary
CDES = "#d9541a"                         # design-point marker line
NEWMAN_S = 0.5
GOLDEN = (1.0 + 5.0 ** 0.5) / 2.0        # overall figure box = golden ratio

# One series per curve: mat key, label, colour, marker
SERIES = [
    ("CT_new", r"$C_T$",                       "#1f77b4", "o"),
    ("CP_new", r"$C_P$",                       "#2ca02c", "s"),
    ("eta",    r"$\eta_{\mathrm{propulsive}}$", "#d62728", "D"),
    ("eps",    r"$\varepsilon$",               "#7f4fbf", "^"),
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
    counts = np.atleast_1d(d["counts"]).astype(float)
    sol = np.atleast_1d(d["solidity"]).astype(float)
    Ndes = int(np.atleast_1d(d["Ndes"])[0])

    fig, ax = plt.subplots(figsize=(8.2, 8.2 / GOLDEN))    # overall box = golden ratio
    fig.subplots_adjust(left=0.10, right=0.965, top=0.86, bottom=0.13)

    # Newman s = 0.5 reference line
    ax.axvline(NEWMAN_S, color=CNEW, ls="--", lw=1.4, zorder=1)

    for key, lab, col, mk in SERIES:
        y = np.atleast_1d(d[key]).astype(float)
        ax.plot(sol, y, "-", lw=2.0, color=col, zorder=2)
        ax.plot(sol, y, mk, ms=5.5, mfc=col, mec=col, zorder=3, label=lab)

    ax.grid(True, alpha=0.25, lw=0.7)
    ax.set_axisbelow(True)
    dx = 0.5 * (sol[1] - sol[0])
    ax.set_xlim(min(sol[0] - dx, NEWMAN_S - 0.1), sol[-1] + dx)
    ax.set_ylim(0.0, 0.75)
    ax.set_xlabel(r"Rotor solidity,  $s = C/S$", fontsize=14)
    ax.set_ylabel(r"[ - ]", fontsize=14)
    ax.tick_params(labelsize=14)

    ax.annotate(r"Newman $s=0.5$", xy=(NEWMAN_S, 0.03), xytext=(4, 3),
                textcoords="offset points", rotation=90, va="bottom", ha="left",
                fontsize=14, color=CNEW)

    # blade count on a secondary top axis
    axt = ax.secondary_xaxis("top")
    axt.set_xticks(sol)
    axt.set_xticklabels([f"{int(c)}" for c in counts], fontsize=14)
    axt.set_xlabel(r"blade count,  $N$", fontsize=14, labelpad=3)
    axt.tick_params(length=3)

    # legend on the figure, top-right corner
    ax.legend(loc="lower right", fontsize=14, handlelength=1.6,
              frameon=True, facecolor="white", framealpha=0.85, edgecolor="0.8")

    out = args.out or os.path.join(DEFAULT_OUT_DIR, "PropCoeffs_vs_Solidity.png")
    fig.savefig(out, dpi=200)
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"saved {out}\n      {pdf}")


if __name__ == "__main__":
    main()
