#!/usr/bin/env python3
"""Blade-chord Reynolds number of the ducted-fan rotor and stator rows.

Uses the same definition as automated-blading/ab_SmithChart.m:

    Re = rho * V * c / mu

with static density rho and static temperature T mass-averaged at the row inlet,
V the relative velocity for the rotor and the absolute velocity for the stator,
mu from Sutherland's law (ab_dynamic_viscosity.m), and c the blade true chord at
mid-span. rho/V/T come from the CFD inlet cuts (<design>_rotor_chic.mat /
<design>_stator_chic.mat, from make_df_chic.m); the chord comes from the blade
splines (<design>_B1_spline.mat rotor, _B2_spline.mat stator).

Usage:
    python calc_reynolds.py
    python calc_reynolds.py Phi06DUCTEDFAN --dir ... --geom-dir ...
"""

import argparse
import glob
import json
import os
import re

import numpy as np
import scipy.io as sio
from scipy.spatial.distance import pdist

from plot_chic import mass_average        # reuse the validated mass-average

DEFAULT_TS_DIR = "/Data/Engine_Selector/TURBOSTREAM"
DEFAULT_GEOM_DIR = "/Data/Engine_Selector/Geometry"
DEFAULT_DESIGNS = ["Phi06DUCTEDFAN", "Phi075DUCTEDFAN", "Phi09DUCTEDFAN2"]

# row -> (chic suffix, spline tag, velocity field)  [ducted fan: rotor + stator]
ROWS = {"rotor": ("rotor", "B1", "V_rel"),
        "stator": ("stator", "B2", "V")}

# propeller operating-point names (cf. PropCURTIS.m: Vf = [1 20 40])
OP_NAMES = {1: "Static", 2: "Take-Off", 3: "Cruise"}


def mu_air(T):
    """Dynamic viscosity of air [kg/m/s] from Sutherland's law (ab_dynamic_viscosity.m)."""
    T0, C, mu0 = 524.07, 120.0, 1.827e-5
    Tr = T * 9.0 / 5.0                     # K -> Rankine
    return mu0 * ((0.555 * T0 + C) / (0.555 * Tr + C)) * (Tr / T0) ** 1.5


def true_chord(spline_mat):
    """Mid-span true chord [m] = largest LE->TE distance of the section in (x, rt)."""
    b = sio.loadmat(spline_mat, struct_as_record=False, squeeze_me=True)["b"]
    xrrt = np.asarray(b.xrrt, float)                  # (i_section, j_span, [x r rt])
    r = xrrt[:, :, 1]
    r_mid = 0.5 * (r.min() + r.max())
    j = int(np.argmin(np.abs(r.mean(axis=0) - r_mid)))
    pts = np.column_stack([xrrt[:, j, 0], xrrt[:, j, 2]])   # (x, rt) at mid-span
    return pdist(pts).max()


def mission_rpm(geom_dir, design):
    """Design (shaft) RPM from the mission file."""
    base = re.sub(r"restart\d*$", "", design)
    with open(os.path.join(geom_dir, base + "-mission.json")) as fh:
        return float(json.load(fh)["RPM"])


def _reynolds(chic_mat, spline_mat, vfield, scale=1.0):
    """Re = rho*V*c/mu from a chic cut file, blade spline and velocity field.
    scale multiplies V (used to project to a different shaft speed, V ~ rpm)."""
    m = sio.loadmat(chic_mat, struct_as_record=False, squeeze_me=True)
    inlet = m["inlet"]
    omega = float(m["rpm"]) * 2.0 * np.pi / 60.0
    ro = mass_average(inlet, "ro", omega)
    T = mass_average(inlet, "T", omega)
    V = mass_average(inlet, vfield, omega) * scale
    c = true_chord(spline_mat)
    mu = mu_air(T)
    return dict(Re=ro * V * c / mu, ro=ro, V=V, T=T, c=c, mu=mu,
                rpm=float(m["rpm"]) * scale)


def design_rows(ts_dir, geom_dir, design, design_rpm=False):
    """Yield (row label, Re dict) for a design, dispatching on machine type:
    ducted fan -> rotor + stator at the design point; propeller -> rotor at each
    operating point (_O1/_O2/_O3). If design_rpm, scale the velocity (and Re) of
    the ducted fan to the current mission RPM (V ~ shaft speed)."""
    op_files = sorted(glob.glob(os.path.join(ts_dir, design + "_O*_avg_chic.mat")))
    if not glob.glob(os.path.join(ts_dir, design + "_rotor_chic.mat")) and op_files:
        # propeller: single rotor, one point per operating condition (rpm unchanged)
        base = re.sub(r"restart\d*$", "", design)               # spline name has no 'restart'
        spline = os.path.join(geom_dir, base + "_B1_spline.mat")
        for f in op_files:
            op = int(re.search(r"_O(\d+)_avg_chic", os.path.basename(f)).group(1))
            yield "rotor/%s" % OP_NAMES.get(op, "O%d" % op), _reynolds(f, spline, "V_rel")
    else:
        # ducted fan: rotor and stator at the design point
        scale = 1.0
        if design_rpm:                  # project as-run velocities to the mission RPM
            mr = sio.loadmat(os.path.join(ts_dir, design + "_rotor_chic.mat"),
                             struct_as_record=False, squeeze_me=True)
            scale = mission_rpm(geom_dir, design) / float(mr["rpm"])
        for row, (suffix, tag, vfield) in ROWS.items():
            chic = os.path.join(ts_dir, "%s_%s_chic.mat" % (design, suffix))
            spline = os.path.join(geom_dir, "%s_%s_spline.mat" % (design, tag))
            yield row, _reynolds(chic, spline, vfield, scale=scale)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designs", nargs="*", default=DEFAULT_DESIGNS)
    ap.add_argument("--dir", default=DEFAULT_TS_DIR, help="TURBOSTREAM directory")
    ap.add_argument("--geom-dir", default=DEFAULT_GEOM_DIR, help="Geometry directory")
    ap.add_argument("--design-rpm", action="store_true",
                    help="scale ducted-fan velocities to the current mission RPM "
                         "(first-order estimate, V ~ shaft speed; CFD not re-run)")
    args = ap.parse_args()

    hdr = "{:<33}{:<15}{:>9}{:>9}{:>9}{:>9}{:>9}{:>12}"
    print(hdr.format("design", "row", "rpm", "rho", "V", "T", "c[mm]", "Re"))
    print("-" * 108)
    for design in args.designs:
        try:
            rows = list(design_rows(args.dir, args.geom_dir, design, args.design_rpm))
        except Exception as exc:
            print("{:<33} {}".format(design, exc))
            continue
        for label, d in rows:
            print("{:<33}{:<15}{:>9.0f}{:>9.3f}{:>9.1f}{:>9.1f}{:>9.2f}{:>12.0f}".format(
                design, label, d["rpm"], d["ro"], d["V"], d["T"], d["c"] * 1e3, d["Re"]))


if __name__ == "__main__":
    main()
