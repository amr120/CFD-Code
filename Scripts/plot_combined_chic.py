#!/usr/bin/env python3
"""Combined propeller + ducted-fan characteristic on one figure.

Plots flow coefficient (phi) against total-total pressure rise coefficient (Cptt)
for BOTH ducted-fan and propeller designs on a single axis, every point coloured
by total-total fan efficiency (eta_poly) on one shared colour scale. Each design
gets its own marker shape; the propulsor type is shown in the legend label.

The two families store their throttle/operating points differently:
  * ducted fan - many throttle points per design, files <design>_M<mass>_avg_chic.mat
                 (cf. plot_chic.py), plotted in descending-phi order.
  * propeller  - three operating points per design, files <design>_O1/_O2/_O3_avg_chic.mat
                 (cf. plot_prop_chic.py), pressure rise measured across the rotor.

A design's type is auto-detected from which files exist (*_M*_ -> ducted fan,
*_O*_ -> propeller). All quantities reproduce the MATLAB definitions in
TS/post/ts_calc_chic.m, so no MATLAB is needed to plot. Run make_df_chic.m /
make_prop_chic.m first if the *_chic.mat cut files are missing.

Usage:
    python plot_combined_chic.py                         # default DF + prop designs
    python plot_combined_chic.py --df ASMEPAPERDF2 --prop 11052026PropPhi06restart2
    python plot_combined_chic.py DESIGN1 DESIGN2         # auto-detect each design's type
    python plot_combined_chic.py --dir /path/to/TURBOSTREAM --out combined.pdf
"""

import argparse
import glob
import os
import re

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Thesis face, matching plot_prop_chic.py: Nimbus Roman body, STIX maths (the
# class takes the `times` option, so the document's maths is mathptmx). The
# heavier black axes frame is kept.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Nimbus Roman", "Nimbus Roman No9 L", "Times New Roman",
                   "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "font.size": 13,
    "axes.edgecolor": "black",
    "axes.linewidth": 1.5,
})

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_DF_DESIGNS = ["Phi06DUCTEDFANFIXRPM", "Phi075DUCTEDFANFIXRPM",
                      "Phi09DUCTEDFAN2FIXRPM"]
DEFAULT_PROP_DESIGNS = ["02062026PropPhi06HYBRIDrestart2",
                        "02062026PropPhi06n05restart2",
                        "11052026PropPhi06restart2"]

# Marker shapes cycled across ALL designs (ducted fan and propeller pooled together)
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h", "p"]

# Operating-point names for propeller designs (cf. PropCURTIS.m: Vf = [1 20 40])
OP_NAMES = {1: "Static", 2: "Take-Off", 3: "Cruise"}

# Explicit legend labels for known propeller designs (blade tip-loading family)
PROP_LABELS = {
    "02062026PropPhi06HYBRIDrestart2": "Hybrid",
    "02062026PropPhi06n05restart2": "Vortex",
    "11052026PropPhi06restart2": "Prandtl Tip",
}


# ----------------------------------------------------------------------------
# Legend labels
# ----------------------------------------------------------------------------
def df_label(design):
    """Ducted-fan legend label, e.g. 'Phi075DUCTEDFAN' -> 'DF, $\\phi_m = 0.75$'."""
    m = re.match(r"Phi(\d+)", design)
    if m:
        digits = m.group(1).lstrip("0") or "0"
        return r"DF, $\phi_m = {:g}$".format(float("0." + digits))
    return "DF, " + design


def prop_label(design):
    """Propeller legend label - explicit override if known, else derived."""
    if design in PROP_LABELS:
        base = PROP_LABELS[design]
    else:
        base = re.sub(r"^\d+Prop", "", design)      # drop leading date + 'Prop'
        base = re.sub(r"restart\d*$", "", base)     # drop trailing 'restart'/'restart2'
        base = base or design
    return "Prop, " + base


# ----------------------------------------------------------------------------
# CFD cut averaging (identical to plot_chic.py / plot_prop_chic.py)
# ----------------------------------------------------------------------------
def _cell_corners(v):
    """Average the four nodal corners of every structured cell -> cell-centred."""
    return 0.25 * (v[:-1, :-1] + v[1:, :-1] + v[:-1, 1:] + v[1:, 1:])


def _cell_areas(x, r, rt):
    """Axial, radial and tangential face-area components of each cell
    (diagonal cross-product, as in ts_mass_average.m / ts_area_average.m)."""
    t = rt / r
    t_av = 0.25 * (t[1:, 1:] + t[:-1, :-1] + t[:-1, 1:] + t[1:, :-1])

    x1 = x[1:, 1:] - x[:-1, :-1]
    x2 = x[:-1, 1:] - x[1:, :-1]
    r1 = r[1:, 1:] - r[:-1, :-1]
    r2 = r[:-1, 1:] - r[1:, :-1]
    rt1 = r[1:, 1:] * (t[1:, 1:] - t_av) - r[:-1, :-1] * (t[:-1, :-1] - t_av)
    rt2 = r[:-1, 1:] * (t[:-1, 1:] - t_av) - r[1:, :-1] * (t[1:, :-1] - t_av)

    Ax = 0.5 * (r1 * rt2 - r2 * rt1)
    Ar = 0.5 * (x2 * rt1 - x1 * rt2)
    At = 0.5 * (x1 * r2 - x2 * r1)
    return Ax, Ar, At


def _mass_flux(cut, Ax, Ar, At, omega):
    """Cell mass flux dM = Ax*rovx + Ar*rovr + At*rovt (ts_mass_average.m)."""
    rovx = _cell_corners(np.asarray(cut.rovx, float))
    rovr = _cell_corners(np.asarray(cut.rovr, float))
    ro = np.asarray(cut.ro, float)
    r = np.asarray(cut.r, float)
    rovt = np.asarray(cut.rorvt, float) / r - r * ro * omega
    rovt = _cell_corners(rovt)
    return Ax * rovx + Ar * rovr + At * rovt


def mass_average(cut, varname, omega):
    """Whole-plane (dim=3) mass average of a stored variable."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    dM = _mass_flux(cut, Ax, Ar, At, omega)
    prop = _cell_corners(np.asarray(getattr(cut, varname), float))
    return np.sum(prop * dM) / np.sum(dM)


def area_average(cut, varname):
    """Whole-plane (dim=3) area average of a stored variable."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    mag_dA = np.sqrt(Ax ** 2 + Ar ** 2 + At ** 2)
    prop = _cell_corners(np.asarray(getattr(cut, varname), float))
    return np.sum(prop * mag_dA) / np.sum(mag_dA)


def calc_point(chic_mat):
    """Return (phi, Cptt, eta_poly) for one *_chic.mat cut file (cf. ts_calc_chic.m).

    Works for both families: the ducted-fan cuts span inlet duct -> outlet duct and
    the propeller cuts span the rotor leading edge -> trailing edge, but the
    definitions of phi, Cptt and eta_poly are identical in each case."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    inlet, outlet, rpm = m["inlet"], m["outlet"], float(m["rpm"])
    omega = rpm * 2.0 * np.pi / 60.0
    ga = float(inlet.av.ga)

    # Midspan blade speed from hub/casing radii at the inlet plane
    r_in = np.asarray(inlet.r, float)
    U_mid = 0.5 * (r_in.min() + r_in.max()) * omega

    Vx_in = area_average(inlet, "Vx")
    ro_in = area_average(inlet, "ro")
    Po_in = mass_average(inlet, "Po", omega)
    Po_out = mass_average(outlet, "Po", omega)
    To_in = mass_average(inlet, "To", omega)
    To_out = mass_average(outlet, "To", omega)

    phi = Vx_in / U_mid
    Cptt = (Po_out - Po_in) / (0.5 * ro_in * U_mid ** 2)
    if To_out > To_in:
        eta_poly = ((ga - 1.0) / ga) * (np.log(Po_out / Po_in) / np.log(To_out / To_in))
    else:
        eta_poly = (ga / (ga - 1.0)) * (np.log(To_out / To_in) / np.log(Po_out / Po_in))
    return phi, Cptt, eta_poly


# ----------------------------------------------------------------------------
# Per-design readers
# ----------------------------------------------------------------------------
def read_df(ts_dir, design):
    """Read every throttle point of a ducted-fan design; phi, Cptt, eta sorted by
    descending phi (as ts_plot_chic sorts)."""
    files = sorted(glob.glob(os.path.join(ts_dir, design + "_M*_avg_chic.mat")))
    phi, cptt, eta = [], [], []
    for f in files:
        try:
            p, c, e = calc_point(f)
        except Exception as exc:
            print("  skipping {}: {}".format(os.path.basename(f), exc))
            continue
        phi.append(p); cptt.append(c); eta.append(e)
    phi, cptt, eta = np.array(phi), np.array(cptt), np.array(eta)
    order = np.argsort(phi)[::-1]
    return phi[order], cptt[order], eta[order]


def read_prop(ts_dir, design):
    """Read the _O1/_O2/_O3 operating points of a propeller design; phi, Cptt, eta
    in operating-point order."""
    ops, phi, cptt, eta = [], [], [], []
    for f in sorted(glob.glob(os.path.join(ts_dir, design + "_O*_avg_chic.mat"))):
        mo = re.search(r"_O(\d+)_avg_chic", os.path.basename(f))
        if mo is None:
            continue
        try:
            p, c, e = calc_point(f)
        except Exception as exc:
            print("  skipping {}: {}".format(os.path.basename(f), exc))
            continue
        ops.append(int(mo.group(1))); phi.append(p); cptt.append(c); eta.append(e)
    ops = np.array(ops); order = np.argsort(ops)
    return (np.array(phi)[order], np.array(cptt)[order], np.array(eta)[order])


def detect_kind(ts_dir, design):
    """Auto-detect a design's family from which cut files exist."""
    if glob.glob(os.path.join(ts_dir, design + "_M*_avg_chic.mat")):
        return "df"
    if glob.glob(os.path.join(ts_dir, design + "_O*_avg_chic.mat")):
        return "prop"
    return None


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*",
                    help="designs to auto-detect (ducted fan or propeller). If "
                         "omitted, the default DF + propeller sets are used.")
    ap.add_argument("--df", nargs="*", default=None,
                    help="ducted-fan designs (default: %s)" % " ".join(DEFAULT_DF_DESIGNS))
    ap.add_argument("--prop", nargs="*", default=None,
                    help="propeller designs (default: %s)" % " ".join(DEFAULT_PROP_DESIGNS))
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--out", default="combined_chic.pdf", help="output figure path")
    ap.add_argument("--no-show", action="store_true", help="save only, do not display")
    args = ap.parse_args()

    # Build the (design, kind) work list. Positional args auto-detect; otherwise the
    # --df / --prop lists (falling back to the defaults when nothing is given).
    work = []
    if args.designs:
        for d in args.designs:
            kind = detect_kind(args.dir, d)
            if kind is None:
                print("  no cut files for %s in %s - skipping" % (d, args.dir))
                continue
            work.append((d, kind))
    if args.df is not None:
        work += [(d, "df") for d in args.df]
    if args.prop is not None:
        work += [(d, "prop") for d in args.prop]
    if not work:                                   # nothing requested -> defaults
        work = ([(d, "df") for d in DEFAULT_DF_DESIGNS]
                + [(d, "prop") for d in DEFAULT_PROP_DESIGNS])

    # Read every design first so the efficiency colour scale is shared across all
    data = []
    for design, kind in work:
        phi, cptt, eta = (read_df(args.dir, design) if kind == "df"
                          else read_prop(args.dir, design))
        if len(phi) == 0:
            print("  no points for %s - skipping" % design)
            continue
        label = df_label(design) if kind == "df" else prop_label(design)
        data.append((design, kind, label, phi, cptt, eta))
        print("{} [{}]: {} points".format(design, kind, len(phi)))
        for p, c, e in zip(phi, cptt, eta):
            print("    phi = {:.4f}   Cptt = {:.4f}   eta = {:.4f}".format(p, c, e))

    if not data:
        raise SystemExit("No designs produced any points - check --dir / design names.")

    all_eta = np.concatenate([d[5] for d in data])
    norm = plt.Normalize(vmin=all_eta.min(), vmax=all_eta.max())
    cmap = plt.get_cmap("viridis")

    golden = (1.0 + 5.0 ** 0.5) / 2.0
    # thesis Figure 5.12, set at 0.85\linewidth (372 pt): drawn 7.14 in (514 pt)
    # so the 13 pt lettering reads at 9.4 pt on the page
    width = 7.14
    fig, ax = plt.subplots(figsize=(width, width / golden))

    handles = []
    for i, (design, kind, label, phi, cptt, eta) in enumerate(data):
        marker = MARKERS[i % len(MARKERS)]
        # neutral grey line connecting the points so fill colour = efficiency only
        ax.plot(phi, cptt, "-", color="0.7", lw=1.0, zorder=1)
        ax.scatter(phi, cptt, c=eta, cmap=cmap, norm=norm, marker=marker,
                   s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles.append(mlines.Line2D([], [], color="0.5", marker=marker,
                                     linestyle="-", markeredgecolor="k",
                                     markeredgewidth=1.1, markersize=10, label=label))

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Total-Total Fan Efficiency")

    ax.set_xlabel("Flow Coefficient")
    ax.set_ylabel("Total-Total Pressure Rise Coefficient")
    ax.grid(True)
    #ax.set_box_aspect(1.0 / golden)        # data box follows the golden ratio (w:h = phi:1)

    # Widen the limits beyond the data so the in-axes legend has clear space
    allphi = np.concatenate([d[3] for d in data])
    allcptt = np.concatenate([d[4] for d in data])
    xpad = 0.12 * (allphi.max() - allphi.min())
    ypad = 0.12 * (allcptt.max() - allcptt.min())
    ax.set_xlim(allphi.min() - 2.6 * xpad, allphi.max() + xpad)
    ax.set_ylim(allcptt.min() - 2.4 * ypad, allcptt.max() + ypad)

    # solid white, not matplotlib's default 0.8 frame alpha, so the grid and
    # the characteristics do not show through the box
    leg = ax.legend(handles=handles, frameon=True, loc="lower left",
                    title="Design", facecolor="white", framealpha=1.0)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("saved", args.out)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
