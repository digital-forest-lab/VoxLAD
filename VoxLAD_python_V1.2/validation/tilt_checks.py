"""Instrument not level: checks of the axis found from the scan lines and of the corrected pulse count.

Simulated scans only, so that anyone can run it:  python validation/tilt_checks.py
(The same checks on real RIEGL VZ-400 scans are described in validation/README.md.)

1. A panoramic scanner in a closed room, tilted by a known angle: the head axis is found from the
   returns in file order.  Lines as planes through the axis, and as cones with an offset beam origin.
2. A box of leaves of known density (LAD = 1 m2/m3) at three elevations, scanned by a scanner leaning
   towards it or away from it: LAD with the axis taken as vertical, and with the axis found.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests"))
import simulate                                                    # noqa: E402
from voxlad import leafangle as la                                # noqa: E402
from voxlad.model import Scan, Settings, compute                  # noqa: E402
from voxlad.orientation import estimate_axis, tilt_effect         # noqa: E402


def angle(a, b):
    return float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b)))))))


def axis_recovery():
    print("1. Head axis found from the returns of a simulated panoramic scan (step 0.1 deg, file in millimetres)")
    print("   lines                         tilt set   tilt found   error     +-      cone found  position error found")
    rng = np.random.default_rng(1)
    step = np.radians(0.1)
    pos = np.array([1.0, 2.0, 0.3])
    worst = 0.0
    for name, kw, perr in (("planes through the axis", {}, (0.0, 0.0)),
                           ("cones of 1 deg, origin 4 mm off", {"cone_deg": 1.0, "origin_offset": 0.004}, (0.0, 0.0)),
                           ("planes, position 10 cm off", {}, (0.08, -0.06))):
        for tilt in (0.0, 0.2, 1.0, 5.0, 30.0, 90.0):
            xyz, axis = simulate.scan_room(pos, simulate.tilt_matrix(tilt, 215.0, 33.0), step, rng,
                                           head_range=(0.0, 360.0) if tilt < 30 else (10.0, 100.0), **kw)
            est = estimate_axis(np.round(xyz, 3), pos + np.array([perr[0], perr[1], 0.0]), step)
            err = angle(est.axis, axis)
            worst = max(worst, err)
            print("   %-30s %6.1f   %9.4f   %7.4f  %6.4f   %8.3f    %5.1f cm %s"
                  % (name, tilt, est.tilt_deg, err, est.se_deg, est.cone_deg, 100 * np.linalg.norm(est.offset),
                     "" if est.ok else "  NOT ACCEPTED: " + est.reason))
    print("   largest error %.4f deg" % worst)
    return worst


def lad_recovery():
    print()
    print("2. Box of leaves of LAD 1.00 (64 voxels of 0.5 m, about 8.5 m from the scanner, step 0.04 deg)")
    print("   elevation  tilt   leaning   N_T error   LAD, axis as vertical   LAD, axis found   tilt found")
    rng = np.random.default_rng(2)
    g = la.leaf_angle_table("spherical")
    st = np.radians(0.04)
    rows = []
    for elev in (20.0, 40.0, 60.0):
        c = 8.5 * np.array([np.cos(np.radians(elev)), 0.0, np.sin(np.radians(elev))])
        lo = (np.round(c * 2) / 2 - 1.0).tolist()
        hi = (np.array(lo) + 2.0).tolist()
        settings = Settings(voxel_size=0.5, bounds=lo + hi, voxels="all", leaf_angle="spherical", jensen_fit="exact")
        for tilt, bearing, word in ((0.0, 90.0, "-"), (1.0, 90.0, "towards"), (1.0, 270.0, "away"), (3.0, 90.0, "towards"), (3.0, 270.0, "away")):
            R = simulate.tilt_matrix(tilt, bearing)
            xyz, _ = simulate.scan_box([0, 0, 0], lo, hi, 1.0, g, st, st, rng, rotation=R)
            cls = np.ones(len(xyz), int)
            out = []
            for axis in ("vertical", "auto"):
                res = compute([Scan(xyz, cls, [0, 0, 0], "spherical", st, axis=axis)], settings)
                out.append((float(np.nanmean(res.lad)), res.scan_axes[0]))
            sign = -1.0 if word == "towards" else 1.0
            nt_err = sign * tilt_effect(tilt, elev) if tilt else 0.0          # leaning towards the box: fewer pulses assumed than fired
            rows.append((elev, tilt, out[0][0], out[1][0]))
            print("   %5.0f     %4.1f   %-8s  %+6.1f %%        %6.3f              %6.3f          %s"
                  % (elev, tilt, word, 100 * nt_err, out[0][0], out[1][0],
                     ("%.2f deg" % out[1][1]["tilt_deg"]) if out[1][1]["source"] == "estimated" else "not found"))
    rows = np.array(rows)
    level = rows[rows[:, 1] == 0]
    tilted = rows[rows[:, 1] > 0]
    ref = {e: v for e, v in zip(level[:, 0], level[:, 3])}
    dev_wrong = max(abs(r[2] / ref[r[0]] - 1) for r in tilted)
    dev_right = max(abs(r[3] / ref[r[0]] - 1) for r in tilted)
    print("   largest change of the mean LAD from the level scan: %.1f %% with the axis as vertical, %.1f %% with the axis found"
          % (100 * dev_wrong, 100 * dev_right))
    return dev_wrong, dev_right


if __name__ == "__main__":
    t = time.time()
    worst = axis_recovery()
    dev_wrong, dev_right = lad_recovery()
    print()
    ok = worst < 0.01 and dev_right < 0.03
    print("%s  (%.0f s)" % ("ALL WITHIN TOLERANCE" if ok else "OUT OF TOLERANCE", time.time() - t))
    sys.exit(0 if ok else 1)
