#!/usr/bin/env python3
"""Ducted-fan post-processing figures from DuctedFanCURTIS.

Figure 1 - the constant-speed characteristic of flow coefficient (phi) against
  total-total pressure rise coefficient (Cptt), coloured by total-total fan
  efficiency (eta_poly). Reads <design>_M<mass>_avg_chic.mat throttle-line cuts.
Figure 2 - the spanwise rotor and stator loss (stagnation pressure loss
  coefficient Yp vs cumulative mass-flow fraction) at the design point. Reads
  <design>_rotor_chic.mat / <design>_stator_chic.mat (run make_df_chic.m first).

All quantities reproduce the MATLAB definitions (TS/post/ts_calc_chic.m,
entrop_loss_prop.m) so no MATLAB is needed to plot. Outputs are <out> and
<out>_loss_span.<ext>.

Usage:
    python plot_chic.py                       # default designs, show + save
    python plot_chic.py Phi06DUCTEDFAN Phi09DUCTEDFAN2
    python plot_chic.py --dir /path/to/TURBOSTREAM --out chic.png
"""

import argparse
import glob
import json
import os
import re

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# Use a 14 pt base font for every text element (labels, ticks, legend, colorbar)
# and a thicker black axes frame.
plt.rcParams.update({"font.size": 14, "axes.edgecolor": "black", "axes.linewidth": 1.5})

# Default location of the TURBOSTREAM solutions (dr.ts in eVTOL.m for this host)
DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_GEOM_DIR = "/Data/Engine_Selector/Geometry"
DEFAULT_DESIGNS = ["Phi06DUCTEDFANFIXRPM", "Phi075DUCTEDFANFIXRPM", "Phi09DUCTEDFAN2FIXRPM"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
ROWS = [("rotor", "-"), ("stator", "--")]      # blade row -> line style (figure 2)
ROW_LOSS_KEY = {"rotor": "RLOSS", "stator": "SLOSS"}   # 2D estimate field per row


def phi_label(design):
    """Legend label from a design name, e.g. 'Phi075DUCTEDFAN' -> '$\\phi_m = 0.75$'."""
    m = re.match(r"Phi(\d+)", design)
    if m:
        digits = m.group(1).lstrip("0") or "0"
        return r"$\phi_m = {:g}$".format(float("0." + digits))
    return design


def _cell_corners(v):
    """Average the four nodal corners of every structured cell -> cell-centred."""
    return 0.25 * (v[:-1, :-1] + v[1:, :-1] + v[:-1, 1:] + v[1:, 1:])


def _cell_areas(x, r, rt):
    """Axial, radial and tangential face-area components of each cell.

    Mirrors the diagonal cross-product area calculation in ts_mass_average.m /
    ts_area_average.m for a structured cut.
    """
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


def _cell_to_node(cell):
    """Spread a radial cell-centred profile back onto the nodal radii
    (matches the re-assignment in ts_mass_average.m / ts_area_average.m)."""
    node = np.empty(len(cell) + 1)
    node[1:-1] = 0.5 * (cell[:-1] + cell[1:])
    node[0] = cell[0]
    node[-1] = cell[-1]
    return node


def radial_profile(cut, varname, omega, weight):
    """Pitchwise (dim=2) average of a variable -> radial nodal profile + per-band
    mass flux. weight='mass' or 'area' (cf. ts_mass/area_average dim=2)."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    dM = _mass_flux(cut, Ax, Ar, At, omega)
    w = dM if weight == "mass" else np.sqrt(Ax ** 2 + Ar ** 2 + At ** 2)
    prop = _cell_corners(np.asarray(getattr(cut, varname), float))
    cell = np.sum(prop * w, axis=1) / np.sum(w, axis=1)
    return _cell_to_node(cell), np.sum(dM, axis=1)


def radial_yp(chic_mat):
    """Spanwise loss profile (cf. entrop_loss_prop.m): cumulative mass-flow
    fraction (0=hub, 1=tip) and Yp at each radial station. Uses relative
    stagnation pressure - for a stationary row (rpm=0) Po_rel = Po automatically."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    inlet, outlet = m["inlet"], m["outlet"]
    omega = float(m["rpm"]) * 2.0 * np.pi / 60.0

    P01, m1 = radial_profile(inlet, "Po_rel", omega, "mass")
    p1, _ = radial_profile(inlet, "P", omega, "area")
    P02, m2 = radial_profile(outlet, "Po_rel", omega, "mass")

    f1 = np.concatenate([[0.0], np.cumsum(m1)]); f1 /= f1[-1]
    f2 = np.concatenate([[0.0], np.cumsum(m2)]); f2 /= f2[-1]
    P01i = np.interp(f2, f1, P01)
    p1i = np.interp(f2, f1, p1)
    return f2, (P01i - P02) / (P01i - p1i)


def radius_massfrac(chic_mat):
    """Map radius -> cumulative mass-flow fraction (0=hub, 1=tip) using the outlet
    plane mass distribution, so a radius-based 2D estimate can share the CFD axis."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    outlet = m["outlet"]
    omega = float(m["rpm"]) * 2.0 * np.pi / 60.0
    x = np.asarray(outlet.x, float)
    r = np.asarray(outlet.r, float)
    rt = np.asarray(outlet.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    band = np.sum(_mass_flux(outlet, Ax, Ar, At, omega), axis=1)
    f = np.concatenate([[0.0], np.cumsum(band)]); f /= f[-1]
    r_node = np.mean(r, axis=1)
    order = np.argsort(r_node)
    return r_node[order], f[order]


def loss_2d_df(geom_dir, design, row):
    """2D cascade-tool loss estimate for a ducted-fan row from the mission file:
    returns (radii, LOSS) using RLOSS/SLOSS over Rh..Rc, or None if absent.
    (The mission stores no radius array, so radii are taken hub->tip uniformly.)"""
    path = os.path.join(geom_dir, design + "-mission.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        mm = json.load(fh)
    key = ROW_LOSS_KEY[row]
    if key not in mm or "Rh" not in mm or "Rc" not in mm:
        return None
    loss = np.asarray(mm[key], float)
    R = np.linspace(float(mm["Rh"]), float(mm["Rc"]), len(loss))
    return R, loss


def calc_point(chic_mat):
    """Return (phi, Cptt, eta_poly) for one *_chic.mat cut file (cf. ts_calc_chic.m)."""
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
    P_dyn = 0.5 * ro_in * U_mid ** 2
    Cptt = (Po_out - Po_in) / P_dyn

    # Total-total polytropic efficiency (ts_calc_chic.m), compressing branch
    if To_out > To_in:
        eta_poly = ((ga - 1.0) / ga) * (np.log(Po_out / Po_in) / np.log(To_out / To_in))
    else:
        eta_poly = (ga / (ga - 1.0)) * (np.log(To_out / To_in) / np.log(Po_out / Po_in))
    return phi, Cptt, eta_poly


def design_chic(ts_dir, design):
    """Read every throttle point for a design; return phi, Cptt, eta sorted by phi."""
    pattern = os.path.join(ts_dir, design + "_M*_avg_chic.mat")
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError("No cut files match " + pattern)

    phi, cptt, eta = [], [], []
    for f in sorted(files):
        try:
            p, c, e = calc_point(f)
        except Exception as exc:  # skip an unreadable / incomplete point
            print("  skipping {}: {}".format(os.path.basename(f), exc))
            continue
        phi.append(p)
        cptt.append(c)
        eta.append(e)

    phi = np.array(phi)
    cptt = np.array(cptt)
    eta = np.array(eta)
    order = np.argsort(phi)[::-1]          # descending phi, as ts_plot_chic sorts
    return phi[order], cptt[order], eta[order]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*", default=DEFAULT_DESIGNS,
                    help="design names (default: %s)" % " ".join(DEFAULT_DESIGNS))
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--geom-dir", default=DEFAULT_GEOM_DIR,
                    help="Geometry directory holding the *-mission.json files")
    ap.add_argument("--out", default="ducted_fan_chic.png", help="output PNG path")
    ap.add_argument("--no-show", action="store_true", help="save only, do not display")
    args = ap.parse_args()

    # Read every design first so the efficiency colour scale is shared
    data = []
    for design in args.designs:
        phi, cptt, eta = design_chic(args.dir, design)
        data.append((design, phi, cptt, eta))
        print("{}: {} points".format(design, len(phi)))
        for p, c, e in zip(phi, cptt, eta):
            print("    phi = {:.4f}   Cptt = {:.4f}   eta = {:.4f}".format(p, c, e))

    all_eta = np.concatenate([d[3] for d in data])
    norm = plt.Normalize(vmin=all_eta.min(), vmax=all_eta.max())
    cmap = plt.get_cmap("viridis")

    import matplotlib.lines as mlines
    golden = (1.0 + 5.0 ** 0.5) / 2.0          # 1.618...
    width = 8.5
    fig, ax = plt.subplots(figsize=(width, width / golden))
    handles = []
    for i, (design, phi, cptt, eta) in enumerate(data):
        marker = MARKERS[i % len(MARKERS)]
        # neutral line connecting the throttle points so colour = efficiency only
        ax.plot(phi, cptt, "-", color="0.7", lw=1.0, zorder=1)
        ax.scatter(phi, cptt, c=eta, cmap=cmap, norm=norm, marker=marker,
                   s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles.append(mlines.Line2D([], [], color="0.5", marker=marker,
                                     linestyle="-", markeredgecolor="k",
                                     markeredgewidth=1.1, markersize=10,
                                     label=phi_label(design)))

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Total-Total Fan Efficiency")

    ax.set_xlabel("Flow Coefficient")
    ax.set_ylabel("Total-Total Pressure Rise Coefficient")
    ax.grid(True)
    leg = ax.legend(handles=handles, frameon=True, loc="lower right")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("saved", args.out)

    # ---- Figure 2: spanwise rotor and stator loss, one panel per design ----
    # Needs <design>_rotor_chic.mat / <design>_stator_chic.mat from make_df_chic.m
    colours = plt.get_cmap("tab10").colors
    row_col = {"rotor": colours[0], "stator": colours[1]}
    n = len(args.designs)
    w2 = 1.0 + 4.0 * n
    fig2, axes = plt.subplots(1, n, sharey=True, figsize=(w2, w2 / golden))
    if n == 1:
        axes = [axes]
    for ax2, design in zip(axes, args.designs):
        for row, _ in ROWS:
            files = glob.glob(os.path.join(args.dir, "%s_%s_chic.mat" % (design, row)))
            if not files:
                print("  no %s cut for %s - run make_df_chic.m" % (row, design))
                continue
            frac, yp = radial_yp(files[0])
            ax2.plot(yp, frac, "-", color=row_col[row], lw=1.6, zorder=3,
                     label=row.capitalize())                       # CFD (solid)
            est = loss_2d_df(args.geom_dir, design, row)
            if est is not None:                                    # 2D estimate (dashed)
                R, loss = est
                r_node, f_node = radius_massfrac(files[0])
                ax2.plot(loss, np.interp(R, r_node, f_node), "--",
                         color=row_col[row], lw=1.6, zorder=2)
            else:
                print("  no %s in %s mission" % (ROW_LOSS_KEY[row], design))
        ax2.set_title(phi_label(design))
        ax2.set_xlim(-0.01, 0.1)            # clip hub/tip end singularities
        ax2.set_ylim(0, 1)
        ax2.grid(True)
    axes[0].set_ylabel("Cumulative Mass-Flow Fraction (hub $\\to$ tip)")
    fig2.supxlabel(r"Stagnation Pressure Loss Coefficient, $Y_p$")
    leg2 = axes[0].legend(frameon=True, loc="center right", title="Row", framealpha=1.0)
    leg2.get_frame().set_edgecolor("black")
    leg2.get_frame().set_facecolor("white")
    leg2.get_frame().set_linewidth(1.3)
    axes[0].add_artist(leg2)
    src_handles = [mlines.Line2D([], [], color="0.3", lw=2, linestyle="-",
                                 label="CFD"),
                   mlines.Line2D([], [], color="0.3", lw=2, linestyle="--",
                                 label="2D estimate")]
    leg2s = axes[-1].legend(handles=src_handles, frameon=True, loc="center right",
                            title="Source", framealpha=1.0)
    leg2s.get_frame().set_edgecolor("black")
    leg2s.get_frame().set_facecolor("white")
    leg2s.get_frame().set_linewidth(1.3)
    fig2.tight_layout()
    base, ext = os.path.splitext(args.out)
    out2 = base + "_loss_span" + ext
    fig2.savefig(out2, dpi=130)
    print("saved", out2)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
