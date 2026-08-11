#!/usr/bin/env python3
"""Downstream sweep of the propeller capture-streamtube jet velocity Vj.

Shows how Vj evolves with axial distance behind the rotor: the slipstream first
contracts (Vj rises toward the ideal jet), then mixes out numerically in the far
field (Vj decays back toward the freestream). Use it to pick a sensible jet
station for eps = Vfs / Vj instead of the mixed-out domain exit.

Reads, per operating point:
  <design>_O#_avg_jetsweep.mat  - constant-x planes from make_prop_jet_sweep.m
  <design>_O#_avg_chic.mat      - rotor mass flow (inlet) and freestream (dom_in)

The streamtube average at each plane reuses streamtube_vx / axial_massflow from
plot_prop_chic.py, so the definition matches the characteristic plots exactly.

Usage:
    python plot_jet_sweep.py                         # default designs
    python plot_jet_sweep.py 02062026PropPhi06HYBRIDrestart2
    python plot_jet_sweep.py --dir /path/to/TURBOSTREAM --out jet_sweep.pdf
"""

import argparse
import glob
import os
import re
import importlib.util

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Reuse the velocity-ratio helpers from plot_prop_chic.py so definitions match
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("ppc", os.path.join(_here, "plot_prop_chic.py"))
ppc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ppc)

plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

OP_NAMES = ppc.OP_NAMES
COND_MARKERS = ppc.COND_MARKERS


def sweep_design(ts_dir, design):
    """Return a list of per-op dicts with the streamtube Vj along the sweep."""
    results = []
    for f in sorted(glob.glob(os.path.join(ts_dir, design + "_O*_avg_chic.mat"))):
        mo = re.search(r"_O(\d+)_avg_chic", os.path.basename(f))
        if mo is None:
            continue
        op = int(mo.group(1))
        sweepf = f.replace("_chic.mat", "_jetsweep.mat")
        if not os.path.exists(sweepf):
            print("  no %s - run make_prop_jet_sweep.m" % os.path.basename(sweepf))
            continue

        c = sio.loadmat(f, struct_as_record=False, squeeze_me=True)
        mdot = ppc.axial_massflow(c["inlet"])
        Vfs = ppc.area_average(c["dom_in"], "Vx")
        Vx_in = ppc.area_average(c["inlet"], "Vx")
        Vx_out = ppc.area_average(c["outlet"], "Vx")

        s = sio.loadmat(sweepf, struct_as_record=False, squeeze_me=True)
        planes = np.atleast_1d(s["planes"])
        xoR = (np.atleast_1d(s["xs"]).astype(float) - float(s["x_te"])) / float(s["Rc"])
        Vj = []
        for pl in planes:
            try:
                Vj.append(ppc.streamtube_vx(pl, mdot))
            except Exception as exc:
                print("    plane skipped: {}".format(exc))
                Vj.append(np.nan)
        results.append(dict(op=op, xoR=xoR, Vj=np.array(Vj),
                            Vfs=Vfs, Vx_in=Vx_in, Vx_out=Vx_out))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*", default=ppc.DEFAULT_DESIGNS,
                    help="design names (default: %s)" % " ".join(ppc.DEFAULT_DESIGNS))
    ap.add_argument("--dir", default=ppc.DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--out", default="jet_sweep.pdf", help="output figure path")
    ap.add_argument("--no-show", action="store_true", help="save only, do not display")
    args = ap.parse_args()

    golden = (1.0 + 5.0 ** 0.5) / 2.0
    width = 8.5
    colours = plt.get_cmap("tab10").colors
    base, ext = os.path.splitext(args.out)

    for design in args.designs:
        data = sweep_design(args.dir, design)
        if not data:
            continue
        label = ppc.design_label(design)
        print("{}:".format(label))

        fig, ax = plt.subplots(figsize=(width, width / golden))
        for d in data:
            col = colours[(d["op"] - 1) % len(colours)]
            ax.plot(d["xoR"], d["Vj"], "-", color=col, marker=COND_MARKERS.get(d["op"], "o"),
                    ms=5, lw=1.6, zorder=3, label=OP_NAMES.get(d["op"], "O%d" % d["op"]))
            # reference levels: rotor-exit (dashed) and freestream (dotted)
            ax.axhline(d["Vx_out"], color=col, ls="--", lw=1.0, alpha=0.6, zorder=1)
            ax.axhline(d["Vfs"], color=col, ls=":", lw=1.0, alpha=0.6, zorder=1)
            # console table so a station can be picked from the numbers
            print("  {}: Vfs={:.2f}  Vx_out={:.2f}".format(
                OP_NAMES.get(d["op"], "O%d" % d["op"]), d["Vfs"], d["Vx_out"]))
            for xo, vj in zip(d["xoR"], d["Vj"]):
                eps = d["Vfs"] / vj if vj else np.nan
                print("      x/Rc = {:6.2f}   Vj = {:7.3f}   eps = {:.4f}".format(xo, vj, eps))

        ax.set_xlabel(r"Axial distance behind rotor TE, $(x - x_{TE}) / R_c$")
        ax.set_ylabel(r"Capture-streamtube jet velocity, $V_j$ (m/s)")
        ax.grid(True)
        leg = ax.legend(frameon=True, title="Condition", loc="best")
        leg.get_frame().set_edgecolor("black")
        leg.get_frame().set_linewidth(1.3)
        ref = [mlines.Line2D([], [], color="0.3", ls="--", label=r"$V_{x,out}$ (rotor exit)"),
               mlines.Line2D([], [], color="0.3", ls=":", label=r"$V_{fs}$ (freestream)")]
        ax.add_artist(leg)
        leg2 = ax.legend(handles=ref, frameon=True, loc="center right", title="Reference")
        leg2.get_frame().set_edgecolor("black")
        leg2.get_frame().set_linewidth(1.3)
        fig.tight_layout()
        out = "{}_{}{}".format(base, ppc.design_label(design).replace(" ", ""), ext)
        fig.savefig(out, dpi=130)
        print("saved", out)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
