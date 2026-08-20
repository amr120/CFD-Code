#!/usr/bin/env python3
"""Plot the propeller characteristic - flow coefficient (phi) against total-total
pressure rise coefficient (Cptt) - with every point coloured by total-total fan
efficiency (eta_poly), in the same style as the ducted-fan plot (plot_chic.py).

Propeller cases differ from the ducted fan:
  * each design has only THREE operating points, ending in _O1/_O2/_O3,
  * there is no duct, so the pressure rise is measured ACROSS THE ROTOR
    (leading-edge face -> trailing-edge face, over the blade span).

The rotor cuts are produced by the MATLAB helper make_prop_chic.m, which writes
<design>_O#_avg_chic.mat holding the inlet/outlet cut planes. This script reads
those and reproduces phi, Cptt and eta_poly with the same definitions as
TS/post/ts_calc_chic.m. Run make_prop_chic.m first if the *_chic.mat are missing.

Usage:
    python plot_prop_chic.py                          # default design
    python plot_prop_chic.py 02062026PropPhi06n05restart2
    python plot_prop_chic.py --dir /path/to/TURBOSTREAM --out prop.pdf
"""

import argparse
import glob
import json
import os
import re

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# Thesis face: Nimbus Roman is the body font, and the class takes the `times`
# option so the maths is mathptmx, which STIX matches. The heavier black axes
# frame is kept. 13 pt is the base; see W_FULL / W_085 for how that lands at
# 9.4 pt on the page.
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

# Thesis \linewidth is 437.46 pt. A figure set at \linewidth is drawn 8.4 in
# (605 pt) wide so its 13 pt lettering reads at 9.4 pt; one set at
# 0.85\linewidth is drawn 7.14 in (514 pt) for the same 9.4 pt.
W_FULL, W_085 = 8.4, 7.14

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_GEOM_DIR = "/Data/Engine_Selector/Geometry"
DEFAULT_DESIGNS = ["02062026PropPhi06HYBRIDrestart2",
                   "02062026PropPhi06n05restart2",
                   "11052026PropPhi06restart2"]

# Operating-point names (cf. leg in PropCURTIS.m: Vf = [1 20 40])
OP_NAMES = {1: "Static", 2: "Take-Off", 3: "Cruise"}

# Marker per operating condition: Static = circle, Take-Off = triangle, Cruise = square
COND_MARKERS = {1: "o", 2: "^", 3: "s"}

# Explicit legend labels for known designs (blade tip-loading family)
LABELS = {
    "02062026PropPhi06HYBRIDrestart2": "Hybrid",
    "02062026PropPhi06n05restart2": "Vortex",
    "11052026PropPhi06restart2": "Prandtl Tip",
}


def design_label(design):
    """Legend label for a design - explicit override if known, else derived, e.g.
    '02062026PropPhi06HYBRIDrestart2' -> 'Hybrid'."""
    if design in LABELS:
        return LABELS[design]
    label = re.sub(r"^\d+Prop", "", design)     # drop leading date + 'Prop'
    label = re.sub(r"restart\d*$", "", label)   # drop trailing 'restart'/'restart2'
    return label or design


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
    """Pitchwise (dim=2) average of a variable -> a radial nodal profile, plus the
    per-band mass flux. weight='mass' or 'area' (cf. ts_mass/area_average dim=2)."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    dM = _mass_flux(cut, Ax, Ar, At, omega)
    w = dM if weight == "mass" else np.sqrt(Ax ** 2 + Ar ** 2 + At ** 2)
    prop = _cell_corners(np.asarray(getattr(cut, varname), float))
    cell = np.sum(prop * w, axis=1) / np.sum(w, axis=1)     # average over pitch (k)
    band_mass = np.sum(dM, axis=1)
    return _cell_to_node(cell), band_mass


def radial_yp(chic_mat):
    """Spanwise rotor loss profile (cf. entrop_loss_prop.m): returns cumulative
    mass-flow fraction (0=hub, 1=tip) and Yp at each radial station."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    inlet, outlet = m["inlet"], m["outlet"]
    omega = float(m["rpm"]) * 2.0 * np.pi / 60.0

    P01, m1 = radial_profile(inlet, "Po_rel", omega, "mass")
    p1, _ = radial_profile(inlet, "P", omega, "area")
    P02, m2 = radial_profile(outlet, "Po_rel", omega, "mass")

    f1 = np.concatenate([[0.0], np.cumsum(m1)]); f1 /= f1[-1]   # inlet mass fraction
    f2 = np.concatenate([[0.0], np.cumsum(m2)]); f2 /= f2[-1]   # outlet mass fraction
    P01i = np.interp(f2, f1, P01)
    p1i = np.interp(f2, f1, p1)
    Yp = (P01i - P02) / (P01i - p1i)
    return f2, Yp


def loss_2d(geom_dir, design):
    """2D cascade-tool estimate of spanwise Yp from the mission file (cf.
    entrop_loss_prop.m): returns the radii R and the LOSS profile, or None if the
    mission file / fields are missing."""
    base = re.sub(r"restart\d*$", "", design)            # strip 'restart'/'restart2'
    path = os.path.join(geom_dir, base + "-mission.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        mm = json.load(fh)
    if "LOSS" not in mm or "R" not in mm:
        return None
    return np.asarray(mm["R"], float), np.asarray(mm["LOSS"], float)


def radius_massfrac(chic_mat):
    """Map radius -> cumulative mass-flow fraction (0=hub, 1=tip) using the outlet
    plane mass distribution, so the 2D estimate can share the CFD mass-fraction axis."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    outlet = m["outlet"]
    omega = float(m["rpm"]) * 2.0 * np.pi / 60.0
    x = np.asarray(outlet.x, float)
    r = np.asarray(outlet.r, float)
    rt = np.asarray(outlet.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    band = np.sum(_mass_flux(outlet, Ax, Ar, At, omega), axis=1)
    f = np.concatenate([[0.0], np.cumsum(band)]); f /= f[-1]
    r_node = np.mean(r, axis=1)                          # nodal radius (mean over pitch)
    order = np.argsort(r_node)                           # ensure increasing for interp
    return r_node[order], f[order]


def axial_massflow(cut):
    """Total axial mass flow through a cut plane, sum(rovx * Ax). Uses only the
    axial momentum component (rovx is the absolute-frame axial flux), so it is
    frame-independent - no rpm/omega needed."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, _, _ = _cell_areas(x, r, rt)
    rovx = _cell_corners(np.asarray(cut.rovx, float))
    return float(np.sum(rovx * Ax))


def streamtube_vx(cut, mdot_target):
    """Area-averaged axial velocity over the capture streamtube at a downstream
    plane. The far-field outlet is mostly undisturbed freestream, so instead of
    averaging the whole plane we accumulate the axial mass flux radially from the
    axis outward until it equals the captured rotor mass flow |mdot_target|, then
    area-average Vx over just that inner (contracted-jet) region. The radial band
    that crosses the target is partially included in proportion to its mass."""
    x = np.asarray(cut.x, float)
    r = np.asarray(cut.r, float)
    rt = np.asarray(cut.rt, float)
    Ax, Ar, At = _cell_areas(x, r, rt)
    rovx = _cell_corners(np.asarray(cut.rovx, float))
    mag_dA = np.sqrt(Ax ** 2 + Ar ** 2 + At ** 2)
    Vx = _cell_corners(np.asarray(cut.Vx, float))

    band_m = np.abs(np.sum(rovx * Ax, axis=1))           # |axial mass flux| per radial band
    band_A = np.sum(mag_dA, axis=1)                      # face area per band
    band_VxA = np.sum(Vx * mag_dA, axis=1)               # Vx . area per band
    r_cell = 0.5 * (r[:-1, :] + r[1:, :]).mean(axis=1)   # radial centre of each band
    order = np.argsort(r_cell)                           # accumulate hub -> tip
    band_m, band_A, band_VxA = band_m[order], band_A[order], band_VxA[order]

    target = abs(mdot_target)
    cum = np.cumsum(band_m)
    if target >= cum[-1]:                                # streamtube fills the whole plane
        return float(np.sum(band_VxA) / np.sum(band_A))
    n = int(np.searchsorted(cum, target))                # first band where cum >= target
    prev = cum[n - 1] if n > 0 else 0.0
    frac = (target - prev) / band_m[n] if band_m[n] > 0 else 0.0   # partial crossing band
    A_inc = band_A[:n].sum() + frac * band_A[n]
    VxA_inc = band_VxA[:n].sum() + frac * band_VxA[n]
    return float(VxA_inc / A_inc)


def jet_velocity(chic_mat, mdot):
    """Developed-jet velocity Vj = the PEAK capture-streamtube Vx over the
    downstream sweep (<file>_jetsweep.mat from make_prop_jet_sweep.m). The
    streamtube contracts and accelerates behind the rotor, peaks at the developed
    jet (~0.8 R_c), then decays as the slipstream mixes out; the peak is the
    physically meaningful jet, independent of any hand-picked station."""
    sweepf = chic_mat.replace("_chic.mat", "_jetsweep.mat")
    if not os.path.exists(sweepf):
        raise FileNotFoundError("%s missing - run make_prop_jet_sweep.m"
                                % os.path.basename(sweepf))
    s = sio.loadmat(sweepf, struct_as_record=False, squeeze_me=True)
    Vj = []
    for pl in np.atleast_1d(s["planes"]):
        try:
            v = streamtube_vx(pl, mdot)
        except Exception:
            continue
        if np.isfinite(v):
            Vj.append(v)
    if not Vj:
        raise ValueError("no valid sweep planes in %s" % os.path.basename(sweepf))
    return float(np.max(Vj))


def calc_point(chic_mat):
    """Return (phi, Cptt, eta_poly, eps, Yp, J, Chi, sigma, Vx_ratio) for one
    *_chic.mat cut file (cf. ts_calc_chic.m / entrop_loss_prop.m).

    Four area-averaged axial velocities define the velocity ratios:
        Vfs    - entering the full domain  (freestream),       from dom_in
        Vj     - developed jet: PEAK capture-streamtube Vx over the downstream
                 sweep (<file>_jetsweep.mat), i.e. where the slipstream is fully
                 contracted but not yet mixed out (~0.8 R_c behind the rotor)
        Vx_in  - entering the blade row    (rotor LE),          from inlet
        Vx_out - exiting the blade row     (rotor TE),          from outlet
    The disk velocity is the mean across the rotor, Vx_avg = 0.5*(Vx_in + Vx_out).
    With U_mid = omega * r_mid (meanline blade speed, the phi denominator):
        phi      = Vx_in  / U_mid          Vx_ratio = Vx_in / Vx_out
        J        = Vfs    / U_mid          Chi      = Vfs / Vx_avg
        eps      = Vfs    / Vj             sigma    = Vx_avg / Vj  (= eps / Chi)
    Yp is the rotor stagnation pressure loss coefficient based on relative
    stagnation pressures: (Po_rel_in - Po_rel_out) / (Po_rel_in - P_in)."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    if "dom_in" not in m:
        raise KeyError("%s has no dom_in - rerun make_prop_chic.m to add the "
                       "full-domain inlet cut" % os.path.basename(chic_mat))
    inlet, outlet, rpm = m["inlet"], m["outlet"], float(m["rpm"])
    dom_in = m["dom_in"]
    omega = rpm * 2.0 * np.pi / 60.0
    ga = float(inlet.av.ga)

    r_in = np.asarray(inlet.r, float)
    r_mid = 0.5 * (r_in.min() + r_in.max())           # meanline radius (as used for phi)
    U_mid = r_mid * omega                             # meanline blade speed (denominator of phi and J)

    Vfs = area_average(dom_in, "Vx")                  # full-domain inlet (freestream)
    Vx_in = area_average(inlet, "Vx")                 # blade-row inlet (rotor LE)
    Vx_out = area_average(outlet, "Vx")               # blade-row outlet (rotor TE)
    mdot_rotor = axial_massflow(inlet)                # captured mass flow through the rotor
    Vj = jet_velocity(chic_mat, mdot_rotor)           # developed jet (peak of the downstream sweep)
    ro_in = area_average(inlet, "ro")
    Po_in = mass_average(inlet, "Po", omega)
    Po_out = mass_average(outlet, "Po", omega)
    To_in = mass_average(inlet, "To", omega)
    To_out = mass_average(outlet, "To", omega)
    Po_rel_in = mass_average(inlet, "Po_rel", omega)
    Po_rel_out = mass_average(outlet, "Po_rel", omega)
    P_in = area_average(inlet, "P")

    Vx_avg = 0.5 * (Vx_in + Vx_out)                   # disk velocity (mean across the rotor)
    phi = Vx_in / U_mid                               # flow coefficient on the blade-row inlet velocity
    Cptt = (Po_out - Po_in) / (0.5 * ro_in * U_mid ** 2)
    J = Vfs / U_mid                                   # advance ratio (same denominator as phi)
    Chi = Vfs / Vx_avg                                # freestream / disk velocity
    eps = Vfs / Vj                                    # freestream / jet (propulsor velocity ratio)
    sigma = Vx_avg / Vj                               # disk velocity / jet (= eps / Chi)
    Vx_ratio = Vx_in / Vx_out                         # axial velocity ratio across the rotor
    Yp = (Po_rel_in - Po_rel_out) / (Po_rel_in - P_in)   # rotor loss coefficient
    if To_out > To_in:
        eta_poly = ((ga - 1.0) / ga) * (np.log(Po_out / Po_in) / np.log(To_out / To_in))
    else:
        eta_poly = (ga / (ga - 1.0)) * (np.log(To_out / To_in) / np.log(Po_out / Po_in))
    return phi, Cptt, eta_poly, eps, Yp, J, Chi, sigma, Vx_ratio


def design_chic(ts_dir, design):
    """Read the _O1/_O2/_O3 operating points; return op numbers, phi, Cptt, eta,
    propulsor efficiency ratio eps = V/V_j, Yp, advance ratio J = V/U_mid, inlet
    velocity ratio Chi = V/Vx_in, sigma = Vx_in/V_j and across-rotor ratio
    Vx_ratio = Vx_in/Vx_out."""
    ops, phi, cptt, eta, eps, yp, J, chi, sig, vxr = [], [], [], [], [], [], [], [], [], []
    for f in sorted(glob.glob(os.path.join(ts_dir, design + "_O*_avg_chic.mat"))):
        mo = re.search(r"_O(\d+)_avg_chic", os.path.basename(f))
        if mo is None:
            continue
        try:
            p, c, e, ep, y, j, ch, sg, vx_ratio = calc_point(f)
        except Exception as exc:
            print("  skipping {}: {}".format(os.path.basename(f), exc))
            continue
        ops.append(int(mo.group(1)))
        phi.append(p)
        cptt.append(c)
        eta.append(e)
        eps.append(ep)                               # V_fs / V_j (freestream / jet)
        yp.append(y)
        J.append(j)                                  # V_fs / U_mid (same denominator as phi)
        chi.append(ch)                               # V_fs / Vx_in (freestream / blade-row inlet)
        sig.append(sg)                               # Vx_in / V_j (blade-row inlet / jet)
        vxr.append(vx_ratio)                         # Vx_in / Vx_out (across the rotor)
    if not ops:
        raise FileNotFoundError(
            "No %s_O*_avg_chic.mat in %s - run make_prop_chic.m first" % (design, ts_dir))
    order = np.argsort(ops)                         # plot in operating-point order
    ops = np.array(ops)[order]
    return (ops, np.array(phi)[order], np.array(cptt)[order],
            np.array(eta)[order], np.array(eps)[order], np.array(yp)[order],
            np.array(J)[order], np.array(chi)[order], np.array(sig)[order],
            np.array(vxr)[order])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*", default=DEFAULT_DESIGNS,
                    help="design names (default: %s)" % " ".join(DEFAULT_DESIGNS))
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--geom-dir", default=DEFAULT_GEOM_DIR,
                    help="Geometry directory holding the *-mission.json files")
    ap.add_argument("--out", default="prop_chic.pdf", help="output figure path")
    ap.add_argument("--span-cond", type=int, default=2, choices=(1, 2, 3),
                    help="operating point for the spanwise Yp profile (1=Static, 2=Take-Off, 3=Cruise)")
    ap.add_argument("--no-show", action="store_true", help="save only, do not display")
    args = ap.parse_args()

    data = []
    for design in args.designs:
        ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio = design_chic(args.dir, design)
        data.append((design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio))
        print("{}:".format(design))
        for o, p, c, e, ep, y, j, ch, sg, vr in zip(ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio):
            print("    O{} {:9s} phi = {:.4f}   Cptt = {:.4f}   eta = {:.4f}   "
                  "eps = {:.4f}   Yp = {:.4f}   J = {:.4f}   Chi = {:.4f}   "
                  "sigma = {:.4f}   Vx_in/Vx_out = {:.4f}".format(
                      o, OP_NAMES.get(o, ""), p, c, e, ep, y, j, ch, sg, vr))

    all_eta = np.concatenate([d[4] for d in data])
    norm = plt.Normalize(vmin=all_eta.min(), vmax=all_eta.max())
    cmap = plt.get_cmap("viridis")

    import matplotlib.lines as mlines

    # Shared legend handles mapping marker shape -> operating condition
    cond_handles = [mlines.Line2D([], [], color="0.4", marker=COND_MARKERS[o],
                                  linestyle="None", markeredgecolor="k", markersize=9,
                                  label=OP_NAMES[o]) for o in sorted(COND_MARKERS)]

    golden = (1.0 + 5.0 ** 0.5) / 2.0
    width = 8.5
    fig, ax = plt.subplots(figsize=(width, width / golden))

    colours = plt.get_cmap("tab10").colors
    handles = []
    for i, (design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio) in enumerate(data):
        col = colours[i % len(colours)]                      # colour = design
        ax.plot(phi, cptt, "-", color=col, lw=1.3, zorder=1)
        for o, p, c, e in zip(ops, phi, cptt, eta):
            ax.scatter(p, c, c=[e], cmap=cmap, norm=norm,    # fill = efficiency
                       marker=COND_MARKERS.get(o, "o"),      # shape = condition
                       s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles.append(mlines.Line2D([], [], color=col, lw=2.5,
                                     label=design_label(design)))

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Total-Total Fan Efficiency")

    ax.set_xlabel("Flow Coefficient")
    ax.set_ylabel("Total-Total Pressure Rise Coefficient")
    ax.grid(True)
    ax.margins(x=0.13, y=0.08)        # headroom so point annotations fit
    leg = ax.legend(handles=handles, frameon=True, loc="lower left", title="Design")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.3)
    ax.add_artist(leg)
    legc = ax.legend(handles=cond_handles, frameon=True, loc="upper right", title="Condition")
    legc.get_frame().set_edgecolor("black")
    legc.get_frame().set_linewidth(1.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("saved", args.out)

    # ---- Figure 2: epsilon (V_in / V_exit) vs total-total fan efficiency ----
    # thesis Figure 5.13, set at 0.85\linewidth
    fig2, ax2 = plt.subplots(figsize=(W_085, W_085 / golden))
    handles2 = []
    for i, (design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio) in enumerate(data):
        col = colours[i % len(colours)]                      # colour = design
        order = np.argsort(eps)                              # tidy connecting line
        ax2.plot(eps[order], eta[order], "-", color=col, lw=1.3, zorder=1)
        for o, ep, e in zip(ops, eps, eta):
            ax2.scatter(ep, e, color=col, marker=COND_MARKERS.get(o, "o"),
                        s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles2.append(mlines.Line2D([], [], color=col, lw=2.5,
                                      label=design_label(design)))

    # varepsilon, as the thesis body and nomenclature set it, not \epsilon, and
    # italic subscripts to match the caption sitting under it on the same page
    ax2.set_xlabel(r"$\varepsilon = V_{fs} / V_j$")
    ax2.set_ylabel("Total-Total Fan Efficiency")
    ax2.grid(True)
    ax2.margins(x=0.13, y=0.10)
    leg2 = ax2.legend(handles=handles2, frameon=True, loc="upper left", title="Design")
    leg2.get_frame().set_edgecolor("black")
    leg2.get_frame().set_linewidth(1.3)
    ax2.add_artist(leg2)
    legc2 = ax2.legend(handles=cond_handles, frameon=True, loc="lower right", title="Condition")
    legc2.get_frame().set_edgecolor("black")
    legc2.get_frame().set_linewidth(1.3)
    fig2.tight_layout()
    base, ext = os.path.splitext(args.out)
    out2 = base + "_eps_eta" + ext
    fig2.savefig(out2, dpi=130)
    print("saved", out2)

    # ---- Figure 3: rotor stagnation pressure loss coefficient Yp vs flow coeff ----
    fig3, ax3 = plt.subplots(figsize=(width, width / golden))
    handles3 = []
    for i, (design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio) in enumerate(data):
        col = colours[i % len(colours)]                      # colour = design
        order = np.argsort(phi)                              # tidy connecting line
        ax3.plot(phi[order], yp[order], "-", color=col, lw=1.3, zorder=1)
        for o, p, y in zip(ops, phi, yp):
            ax3.scatter(p, y, color=col, marker=COND_MARKERS.get(o, "o"),
                        s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles3.append(mlines.Line2D([], [], color=col, lw=2.5,
                                      label=design_label(design)))

    ax3.set_xlabel("Flow Coefficient")
    ax3.set_ylabel(r"Stagnation Pressure Loss Coefficient, $Y_p$")
    ax3.grid(True)
    ax3.margins(x=0.13, y=0.10)
    leg3 = ax3.legend(handles=handles3, frameon=True, loc="upper right", title="Design")
    leg3.get_frame().set_edgecolor("black")
    leg3.get_frame().set_linewidth(1.3)
    ax3.add_artist(leg3)
    legc3 = ax3.legend(handles=cond_handles, frameon=True, loc="lower left", title="Condition")
    legc3.get_frame().set_edgecolor("black")
    legc3.get_frame().set_linewidth(1.3)
    fig3.tight_layout()
    out3 = base + "_yp" + ext
    fig3.savefig(out3, dpi=130)
    print("saved", out3)

    # ---- Figure 4: spanwise rotor loss profile Yp(span) at one condition ----
    cond = args.span_cond
    # thesis Figure 5.11, set at \linewidth
    fig4, ax4 = plt.subplots(figsize=(W_FULL, W_FULL / golden))
    design_handles4 = []
    for i, design in enumerate(args.designs):
        col = colours[i % len(colours)]
        files = glob.glob(os.path.join(args.dir,
                                       "%s_O%d_avg_chic.mat" % (design, cond)))
        cfd_file = files[0] if files else None
        if cfd_file:
            frac, yp = radial_yp(cfd_file)
            ax4.plot(yp, frac, "-", color=col, lw=1.6, zorder=2)     # CFD (solid)
        else:
            print("  no O%d cut for %s - skipping CFD span profile" % (cond, design))
        est = loss_2d(args.geom_dir, design)
        if est is not None:                                          # 2D estimate (dashed)
            R, loss = est
            if cfd_file is not None:        # map radius -> mass fraction (shared axis)
                r_node, f_node = radius_massfrac(cfd_file)
                yloss = np.interp(R, r_node, f_node)
            else:                           # no CFD to map with -> geometric span
                yloss = (R - R.min()) / (R.max() - R.min())
            ax4.plot(loss, yloss, "--", color=col, lw=1.6, zorder=2)
        else:
            print("  no 2D LOSS in mission for %s" % design)
        design_handles4.append(mlines.Line2D([], [], color=col, lw=2.5,
                                             label=design_label(design)))

    # Second legend distinguishing the two sources (solid CFD vs dashed 2D estimate)
    src_handles = [mlines.Line2D([], [], color="0.3", lw=2, linestyle="-",
                                 label="CFD (TURBOSTREAM)"),
                   mlines.Line2D([], [], color="0.3", lw=2, linestyle="--",
                                 label="2D cascade estimate")]

    ax4.set_xlabel(r"Stagnation Pressure Loss Coefficient, $Y_p$")
    ax4.set_ylabel("Cumulative Mass-Flow Fraction (hub $\\to$ tip)")
    # No on-figure title: the caption of Figure 5.11 names the condition.
    ax4.set_xlim(-0.01, 0.08)            # clip hub/tip end singularities (cf. entrop_loss_prop.m)
    ax4.set_ylim(0, 1)
    ax4.grid(True)
    leg4 = ax4.legend(handles=design_handles4, frameon=True, loc="upper left", title="Design")
    leg4.get_frame().set_edgecolor("black")
    leg4.get_frame().set_linewidth(1.3)
    ax4.add_artist(leg4)
    leg4s = ax4.legend(handles=src_handles, frameon=True, loc="lower left", title="Source")
    leg4s.get_frame().set_edgecolor("black")
    leg4s.get_frame().set_linewidth(1.3)
    fig4.tight_layout()
    out4 = base + "_yp_span" + ext
    fig4.savefig(out4, dpi=130)
    print("saved", out4)

    # ---- Figure 5: advance ratio J vs flow coefficient ----
    fig5, ax5 = plt.subplots(figsize=(width, width / golden))
    handles5 = []
    for i, (design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio) in enumerate(data):
        col = colours[i % len(colours)]                      # colour = design
        order = np.argsort(phi)                              # tidy connecting line
        ax5.plot(phi[order], J[order], "-", color=col, lw=1.3, zorder=1)
        for o, p, j in zip(ops, phi, J):
            ax5.scatter(p, j, color=col, marker=COND_MARKERS.get(o, "o"),
                        s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles5.append(mlines.Line2D([], [], color=col, lw=2.5,
                                      label=design_label(design)))

    ax5.set_xlabel("Flow Coefficient")
    ax5.set_ylabel(r"Advance Ratio, $J = V / U$")
    ax5.grid(True)
    ax5.margins(x=0.13, y=0.10)
    leg5 = ax5.legend(handles=handles5, frameon=True, loc="upper left", title="Design")
    leg5.get_frame().set_edgecolor("black")
    leg5.get_frame().set_linewidth(1.3)
    ax5.add_artist(leg5)
    legc5 = ax5.legend(handles=cond_handles, frameon=True, loc="lower right", title="Condition")
    legc5.get_frame().set_edgecolor("black")
    legc5.get_frame().set_linewidth(1.3)
    fig5.tight_layout()
    out5 = base + "_J_phi" + ext
    fig5.savefig(out5, dpi=130)
    print("saved", out5)

    # ---- Figure 6: Chi (Vfs / Vx_in) vs sigma (Vx_in / Vj), with constant-eps
    #      contours (eps = Vfs/Vj = Chi * sigma) ----
    fig6, ax6 = plt.subplots(figsize=(width, width / golden))

    # Axis limits from the data with a little headroom
    allsig = np.concatenate([d[9] for d in data])
    allchi = np.concatenate([d[8] for d in data])
    sx = 0.10 * (allsig.max() - allsig.min()); cy = 0.10 * (allchi.max() - allchi.min())
    xlim = (max(0.0, allsig.min() - sx), allsig.max() + sx)
    ylim = (max(0.0, allchi.min() - cy), allchi.max() + cy)

    # Constant-eps contours: eps = sigma * chi (hyperbolae chi = eps / sigma)
    sgrid = np.linspace(xlim[0], xlim[1], 250)
    cgrid = np.linspace(ylim[0], ylim[1], 250)
    S, C = np.meshgrid(sgrid, cgrid)
    levels = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    cs = ax6.contour(S, C, S * C, levels=levels, colors="0.6",
                     linewidths=0.8, linestyles="--", zorder=0)
    ax6.clabel(cs, fmt=r"$\epsilon=%.2f$", fontsize=9, inline=True)

    # Actuator-disk line: Vx_avg = 0.5 (Vfs + Vj)  ->  chi = 2 - 1/sigma
    ad_line, = ax6.plot(sgrid, 2.0 - 1.0 / sgrid, "-", color="k", lw=1.8, zorder=2.5,
                        label=r"Actuator disk, $V_{x,\mathrm{avg}}=\frac{1}{2}(V_{\mathrm{fs}}+V_j)$")

    handles6 = []
    for i, (design, ops, phi, cptt, eta, eps, yp, J, Chi, sigma, Vx_ratio) in enumerate(data):
        col = colours[i % len(colours)]                      # colour = design
        order = np.argsort(sigma)                            # tidy connecting line
        ax6.plot(sigma[order], Chi[order], "-", color=col, lw=1.3, zorder=2)
        for o, ch, sg in zip(ops, Chi, sigma):
            ax6.scatter(sg, ch, color=col, marker=COND_MARKERS.get(o, "o"),
                        s=100, edgecolors="k", linewidths=1.1, zorder=3)
        handles6.append(mlines.Line2D([], [], color=col, lw=2.5,
                                      label=design_label(design)))

    ax6.set_xlabel(r"$\sigma = V_{x,\mathrm{avg}} / V_j$")
    ax6.set_ylabel(r"$\chi = V_{\mathrm{fs}} / V_{x,\mathrm{avg}}$")
    ax6.set_xlim(xlim); ax6.set_ylim(ylim)
    ax6.grid(True)
    leg6 = ax6.legend(handles=handles6, frameon=True, loc="upper left", title="Design")
    leg6.get_frame().set_edgecolor("black")
    leg6.get_frame().set_linewidth(1.3)
    ax6.add_artist(leg6)
    legc6 = ax6.legend(handles=cond_handles, frameon=True, loc="lower right", title="Condition")
    legc6.get_frame().set_edgecolor("black")
    legc6.get_frame().set_linewidth(1.3)
    ax6.add_artist(legc6)
    legad6 = ax6.legend(handles=[ad_line], frameon=True, loc="lower center")
    legad6.get_frame().set_edgecolor("black")
    legad6.get_frame().set_linewidth(1.3)
    fig6.tight_layout()
    out6 = base + "_sigma_chi" + ext
    fig6.savefig(out6, dpi=130)
    print("saved", out6)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
