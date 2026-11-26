"""The tilt handling checked on real scans of a RIEGL VZ-400 (not distributed with VoxLAD).

What it needs: the pulses of each raw .rxp file written by the reader `rxp_shots` (RIEGL's RiVLib),
with one line added so that it also writes the beam origin of each pulse (see README.md here):

    <name>_shots.bin, <name>_echoes.bin, <name>_origins.bin

Usage:  python vz400_checks.py FOLDER UPRIGHT_NAME [TILTED_NAME MATRIX_FILE]

UPRIGHT_NAME  a scan with the instrument upright (its own frame is used: the head axis is z exactly)
TILTED_NAME   a scan taken on a tilt mount; MATRIX_FILE holds the 4 x 4 registration matrix that
              places it in the frame of the upright scan

Checks
  A. the head axis found from the returns, after the cloud was turned by known angles, cropped,
     rounded, thinned, shuffled, or given a wrong scanner position;
  B. the angular steps measured from the returns against those of all the pulses of the file;
  C. N_T of the model against the number of pulses the instrument actually fired at 150 voxels,
     with the axis taken as vertical and with the true axis;
  D. the same for the tilt-mount scan placed in the frame of the upright scan.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests"))
import simulate                                                    # noqa: E402
from voxlad import geometry as g                                  # noqa: E402
from voxlad.orientation import estimate_axis                      # noqa: E402

SHOT = np.dtype([("dir", "<f8", 3), ("time", "<f8"), ("n", "u1")])
ECHO = np.dtype([("shot", "<u4"), ("i", "u1"), ("range", "<f8"), ("amp", "<f4"), ("refl", "<f4"), ("dev", "<f4")])


def load(folder, name):
    """Directions and origins of all pulses, and the returns as an exported cloud would hold them."""
    S = np.memmap(os.path.join(folder, name + "_shots.bin"), dtype=SHOT, mode="r")
    E = np.memmap(os.path.join(folder, name + "_echoes.bin"), dtype=ECHO, mode="r")
    d = np.asarray(S["dir"], dtype=np.float64)
    o = np.fromfile(os.path.join(folder, name + "_origins.bin"), dtype="<f4").reshape(-1, 3)
    es, er = np.asarray(E["shot"]), np.asarray(E["range"], dtype=np.float64)
    xyz = np.empty((len(es), 3))
    for a in range(0, len(es), 4_000_000):
        s = es[a:a + 4_000_000]
        xyz[a:a + 4_000_000] = o[s].astype(np.float64) + er[a:a + 4_000_000, None] * d[s]
    return d, o, xyz


def steps_of_pulses(d):
    """Mean angular steps along and between the scan lines, from every pulse of the file."""
    th = np.arccos(np.clip(d[:, 2], -1, 1))
    ph = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    brk = np.flatnonzero(np.abs(np.diff(th)) > np.radians(5)) + 1
    along = float(np.median(np.abs(np.diff(th[brk[3]:brk[4]]))))
    mid = (brk[:-1] + brk[1:]) // 2                      # the mirror facets give unequal spacings: use the mean
    between = float(abs(ph[mid[-1]] - ph[mid[0]]) / (len(mid) - 1))
    return between, along, len(brk) + 1, (ph[mid[0]], ph[mid[-1]])


def count_pulses(d32, o32, R, pos, lo, h):
    """Pulses whose beam crosses each voxel (lower corners lo, side h, in the common frame)."""
    c = lo + h / 2
    out = np.zeros(len(lo), int)
    for i in range(len(lo)):
        v = (c[i] - pos) @ R
        r = np.linalg.norm(v)
        cosa = np.cos(np.arcsin(min(1.0, (h * 0.87 + 0.12) / r)) + 0.002)
        cand = np.flatnonzero(d32 @ (v / r).astype(np.float32) > cosa)
        dd = d32[cand].astype(float) @ R.T
        oo = o32[cand].astype(float) @ R.T + pos
        with np.errstate(divide="ignore", invalid="ignore"):
            ta, tb = (lo[i] - oo) / dd, (lo[i] + h - oo) / dd
        t0, t1 = np.minimum(ta, tb).max(axis=1), np.maximum(ta, tb).min(axis=1)
        out[i] = int(((t1 > t0) & (t1 > 0)).sum())
    return out


def pick_voxels(rng, n, pos, axis, h, limits, inside):
    lo = []
    while len(lo) < n:
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        c = pos + u * rng.uniform(6, 25)
        if g.covered((c - h / 2)[None], h, pos, None, None, limits[0], limits[1], axis=axis)[0] and inside(c - pos):
            lo.append(c - h / 2)
    return np.array(lo)


def nt_report(lo, h, pos, actual, between, along, axes):
    c = lo + h / 2
    ap, _ = g.projected_area(lo, h, pos)
    for name, ax in axes:
        e = g.pulse_count(ap, c, pos, "spherical", between, along, axis=ax) / actual - 1
        print("      %-30s N_T model / pulses fired - 1: mean %+6.2f %%, rms %5.2f %%, worst %+6.1f %%"
              % (name, 100 * e.mean(), 100 * np.sqrt((e ** 2).mean()), 100 * e[np.argmax(np.abs(e))]))


def angle(a, b):
    return float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b)))))))


def main(folder, upright, tilted=None, matrix_file=None):
    rng = np.random.default_rng(3)
    z = np.array([0.0, 0.0, 1.0])
    d, o, xyz = load(folder, upright)
    between, along, nlines, seam = steps_of_pulses(d)
    print("Upright scan %s: %d pulses in %d lines, %d returns; steps %.5f deg between lines (mean), %.5f deg along them"
          % (upright, len(d), nlines, len(xyz), np.degrees(between), np.degrees(along)))
    step = np.radians(0.04)

    print("A. head axis found from the returns (true axis known)")
    rho = np.linalg.norm(xyz, axis=1)
    phi = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])) % 360
    th = np.degrees(np.arccos(xyz[:, 2] / rho))
    R = simulate.tilt_matrix(1.7, 127.0, -71.0)
    true = R @ z
    pos = np.array([5.0, 6.0, 1.5])
    cases = [("all returns", slice(None), (0, 0, 0), 0.001),
             ("sector of 30 deg", (phi > 100) & (phi < 130), (0, 0, 0), 0.001),
             ("sector 10 deg, elevation 10 to 35 deg", (phi > 100) & (phi < 110) & (th < 80) & (th > 55), (0, 0, 0), 0.001),
             ("the same, 8 to 25 m only", (phi > 100) & (phi < 110) & (th < 80) & (th > 55) & (rho > 8) & (rho < 25), (0, 0, 0), 0.001),
             ("sector 4 deg, elevation 20 to 30 deg", (phi > 100) & (phi < 104) & (th < 70) & (th > 60) & (rho > 15) & (rho < 25), (0, 0, 0), 0.001),
             ("box of 6 x 6 m, 10 m away", (np.abs(xyz[:, 0] - 10) < 3) & (np.abs(xyz[:, 1] - 2) < 3), (0, 0, 0), 0.001),
             ("file rounded to centimetres", slice(None), (0, 0, 0), 0.01),
             ("every 4th return kept", slice(None, None, 4), (0, 0, 0), 0.001),
             ("declared position 5 cm off", slice(None), (0.03, -0.04, 0.0), 0.001),
             ("declared position 30 cm off", slice(None), (0.2, -0.2, 0.1), 0.001),
             ("declared position 20 cm too high", slice(None), (0.0, 0.0, 0.2), 0.001),
             ("returns in random order", rng.permutation(len(xyz))[:3_000_000], (0, 0, 0), 0.001)]
    print("   cloud turned by 1.7 deg                 accepted  error (deg)   +-     cone   position error found   steps found (deg)")
    for name, sel, perr, q in cases:
        w = np.round((xyz[sel] @ R.T + pos) / q) * q
        e = estimate_axis(w, pos + np.array(perr), step)
        print("   %-40s %-5s   %8.4f  %6.4f  %5.3f   %6.1f cm            %.5f  %.5f   %s"
              % (name, "yes" if e.ok else "no", angle(e.axis, true), e.se_deg, e.cone_deg, 100 * np.linalg.norm(e.offset),
                 np.degrees(e.step_between), np.degrees(e.step_along), "" if e.ok else e.reason))
    for tilt in (0.0, 0.3, 3.0, 10.0, 45.0, 90.0):
        Rt = simulate.tilt_matrix(tilt, 250.0, 33.0)
        e = estimate_axis(np.round(xyz @ Rt.T + pos, 3), pos, step)
        print("   cloud turned by %4.1f deg: error %.4f deg" % (tilt, angle(e.axis, Rt @ z)))

    print("B. see the last two columns above; from all pulses: %.5f and %.5f deg" % (np.degrees(between), np.degrees(along)))

    full = abs(seam[1] - seam[0]) > np.radians(355) and np.degrees(np.arccos(d[:, 2].max())) < 31 \
        and np.degrees(np.arccos(d[:, 2].min())) > 129
    print("C. N_T against the pulses fired (150 voxels of 0.5 m, 6 to 25 m from the scanner)"
          + ("" if full else ": skipped, the scan is not a full turn from 30 to 130 degrees"))
    d32 = d.astype(np.float32)
    h = 0.5
    limits = (np.radians(30.2), np.radians(129.8))
    s0 = np.degrees(seam[0]) % 360
    for tilt in ((0.0, 1.0, 3.0, 10.0) if full else ()):
        Rt = simulate.tilt_matrix(tilt, 127.0, -71.0)
        ax = Rt @ z

        def away(v, Rt=Rt):                       # not in the few degrees where the scan overlaps itself
            vs = v @ Rt
            a = np.degrees(np.arctan2(vs[1], vs[0])) % 360
            return min(abs(a - s0), 360 - abs(a - s0)) > 8.0
        lo = pick_voxels(rng, 150, pos, ax, h, limits, away)
        t = time.time()
        actual = count_pulses(d32, o, Rt, pos, lo, h)
        print("   frame turned by %.0f deg: %d to %d pulses per voxel (%.0f s)" % (tilt, actual.min(), actual.max(), time.time() - t))
        nt_report(lo, h, pos, actual, between, along, (("axis taken as vertical", None), ("true axis", ax)))
    del d, d32, xyz

    if tilted:
        M = np.loadtxt(matrix_file).reshape(4, 4)
        R, pos = M[:3, :3], M[:3, 3]
        ax = R[:, 2]
        d, o, xyz = load(folder, tilted)
        between, along, nlines, _ = steps_of_pulses(d)
        print("D. Tilt-mount scan %s in the frame of the upright scan: %d pulses, axis %.2f deg from the vertical"
              % (tilted, len(d), angle(ax, z)))
        e = estimate_axis(np.round(xyz @ R.T + pos, 3), pos, step)
        print("   axis found from the returns: accepted %s, %.4f deg from the axis of the registration matrix (+-%.4f)"
              % ("yes" if e.ok else "no", angle(e.axis, ax), e.se_deg))
        th = np.arccos(d[:, 2])
        ph = np.arctan2(d[:, 1], d[:, 0])
        limits = (th.min() + np.radians(0.3), th.max() - np.radians(0.3))
        p0, p1 = ph.min() + np.radians(0.5), ph.max() - np.radians(0.5)

        def inside(v):
            vs = (v[None, :] + (g.VERTEX_OFFSETS - 0.5) * h) @ R
            a = np.arctan2(vs[:, 1], vs[:, 0])
            return bool(np.all((a > p0) & (a < p1)))
        lo = pick_voxels(rng, 150, pos, ax, h, limits, inside)
        actual = count_pulses(d.astype(np.float32), o, R, pos, lo, h)
        nt_report(lo, h, pos, actual, between, along,
                  (("axis taken as vertical", None), ("axis of the registration", ax), ("axis found from the returns", e.axis)))


if __name__ == "__main__":
    if len(sys.argv) not in (3, 5):
        sys.exit(__doc__)
    main(*sys.argv[1:])
