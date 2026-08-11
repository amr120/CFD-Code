#!/usr/bin/env python3
"""Plot the spanwise rotor and stator loss (stagnation pressure loss coefficient
Yp vs cumulative mass-flow fraction) for the ducted-fan design points, in the
same style as the propeller span plot (plot_prop_chic.py figure 4).

The ducted fan has two blade rows. The rotor loss is evaluated with relative
stagnation pressure, the stator with absolute - handled automatically because
make_df_chic.m saves each row's cut with the correct rpm (stator rpm = 0, so
Po_rel = Po). Run make_df_chic.m first to generate <design>_rotor_chic.mat and
<design>_stator_chic.mat.

Usage:
    python plot_df_loss.py
    python plot_df_loss.py Phi06DUCTEDFAN Phi09DUCTEDFAN2 --out df_loss.pdf
"""

import argparse
import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Reuse the validated cut-averaging / loss code from the propeller script
from plot_prop_chic import radial_yp

plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_DESIGNS = ["Phi06DUCTEDFAN", "Phi075DUCTEDFAN", "Phi09DUCTEDFAN2"]
ROWS = [("rotor", "-"), ("stator", "--")]      # row name -> line style


def phi_label(design):
    """Legend label from a design name, e.g. 'Phi075DUCTEDFAN' -> '$\\phi_m = 0.75$'."""
    m = re.search(r"Phi(\d+)", design)
    if m:
        digits = m.group(1).lstrip("0") or "0"
        return r"$\phi_m = {:g}$".format(float("0." + digits))
    return design


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*", default=DEFAULT_DESIGNS,
                    help="design names (default: %s)" % " ".join(DEFAULT_DESIGNS))
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--out", default="df_loss_span.pdf", help="output figure path")
    ap.add_argument("--no-show", action="store_true", help="save only, do not display")
    args = ap.parse_args()

    colours = plt.get_cmap("tab10").colors
    golden = (1.0 + 5.0 ** 0.5) / 2.0
    width = 8.5
    fig, ax = plt.subplots(figsize=(width, width / golden))

    design_handles = []
    for i, design in enumerate(args.designs):
        col = colours[i % len(colours)]
        plotted = False
        for row, style in ROWS:
            files = glob.glob(os.path.join(args.dir, "%s_%s_chic.mat" % (design, row)))
            if not files:
                print("  no %s cut for %s - run make_df_chic.m" % (row, design))
                continue
            frac, yp = radial_yp(files[0])
            ax.plot(yp, frac, style, color=col, lw=1.6, zorder=2)
            plotted = True
        if plotted:
            design_handles.append(mlines.Line2D([], [], color=col, lw=2.5,
                                                label=phi_label(design)))

    # Second legend distinguishing rotor (solid) from stator (dashed)
    row_handles = [mlines.Line2D([], [], color="0.3", lw=2, linestyle=style,
                                 label=row.capitalize()) for row, style in ROWS]

    ax.set_xlabel(r"Stagnation Pressure Loss Coefficient, $Y_p$")
    ax.set_ylabel("Cumulative Mass-Flow Fraction (hub $\\to$ tip)")
    ax.set_xlim(-0.01, 0.1)             # clip hub/tip end singularities
    ax.set_ylim(0, 1)
    ax.grid(True)
    leg = ax.legend(handles=design_handles, frameon=True, loc="upper left", title="Design")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.3)
    ax.add_artist(leg)
    legr = ax.legend(handles=row_handles, frameon=True, loc="lower left", title="Row")
    legr.get_frame().set_edgecolor("black")
    legr.get_frame().set_linewidth(1.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("saved", args.out)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
