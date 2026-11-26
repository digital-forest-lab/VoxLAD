"""Checks on a scan before running the model."""
from __future__ import annotations

import numpy as np

from . import geometry as geo
from . import orientation
from .io import LEAF, WOOD, NOISE


def angles(scan):
    """Bearing (clockwise from +y), zenith angle and range of every return.  For a scanner that
    was not level, the angles are taken about the axis its head turned around, as its limits are."""
    v = scan.xyz - scan.position
    rng = np.sqrt((v ** 2).sum(axis=1))
    zen, az = geo.axis_angles(v, orientation.resolve_axis(scan)["axis"])
    return az, zen, rng


def density_check(scan, window=30):
    """Compare the densest patches of returns with what the declared angular
    step allows.

    The sky is divided into windows ``window`` steps wide.  Where a surface
    returns every pulse (ground, a trunk) a window holds about window**2
    returns; it cannot hold more unless the real step is finer than declared
    or the file holds several returns per pulse.

    Returns (ratio, n_windows): ratio = returns in the densest windows (99.5th
    percentile of occupied windows) / returns the declared step allows.
    """
    az, zen, rng = angles(scan)
    ok = rng > 0
    az, zen = az[ok], zen[ok]
    if az.size < 1000:
        return np.nan, 0
    # unwrap azimuth around the mean direction so a sector is not split at +-180
    mid = np.arctan2(np.mean(np.sin(az)), np.mean(np.cos(az)))
    a = np.mod(az - mid + np.pi, 2 * np.pi) - np.pi
    if scan.pattern == "spherical":
        u, du = a, window * scan.step_az                  # equal steps in azimuth and zenith
        w, dw = zen, window * scan.step_zen
        allowed = float(window) ** 2
        ii = np.floor(u / du).astype(np.int64)
        jj = np.floor(w / dw).astype(np.int64)
    else:
        # the same number of pulses per unit solid angle everywhere: use equal-area cells
        du = window * scan.step_az
        ii = np.floor(a / du).astype(np.int64)
        dw = window * scan.step_zen
        # cells of equal solid angle du*dw: bands of equal width in cos(zenith) at the scan's mean zenith
        s0 = max(np.sin(np.median(zen)), 0.05)
        jj = np.floor((np.cos(zen)) / (dw * s0)).astype(np.int64)
        allowed = float(window) ** 2 / s0
    key = (ii - ii.min()) * (jj.max() - jj.min() + 1) + (jj - jj.min())
    counts = np.bincount(key)
    counts = counts[counts > 0]
    return float(np.percentile(counts, 99.5) / allowed), int(counts.size)


def describe(scan):
    """Text block describing one scan."""
    az, zen, rng = angles(scan)
    L = []
    n = scan.xyz.shape[0]
    L.append(f"scan '{scan.name}': {n} returns used"
             + (f" of {scan.n_in_file} in the file" if getattr(scan, "n_in_file", n) != n else ""))
    L.append(f"  classes        leaf {int(np.sum(scan.cls == LEAF))}, wood {int(np.sum(scan.cls == WOOD))}, "
             f"noise {int(np.sum(scan.cls == NOISE))}")
    L.append(f"  scanner at     ({scan.position[0]:g}, {scan.position[1]:g}, {scan.position[2]:g}); "
             f"pattern {scan.pattern}; step {np.degrees(scan.step_az):.5g} x {np.degrees(scan.step_zen):.5g} deg "
             f"({scan.step_az:.3g} x {scan.step_zen:.3g} rad)")
    if n == 0:
        return "\n".join(L)
    info = orientation.resolve_axis(scan)
    L.append("  head axis      " + info["text"])
    est = info["estimate"]
    if est is not None and est.ok:
        if info["source"] == "estimated" and 0.05 <= info["tilt_deg"] <= 20.0:
            e30, e60 = (100 * orientation.tilt_effect(info["tilt_deg"], e) for e in (30, 60))
            L.append(f"                 taken as vertical, this tilt would put N_T off by up to {e30:.1f}% at 30 deg "
                     f"elevation and {e60:.1f}% at 60 deg; it is accounted for")
        L.append(f"  scan lines     {est.lines:,} lines analysed; cones of {est.cone_deg:.2f} deg; they meet "
                 f"{1000 * float(np.linalg.norm(est.offset)):.0f} mm from the declared position; "
                 f"returns lie {est.rms_mm:.1f} mm from their line")
        sb = f"{np.degrees(est.step_between):.5g}" if np.isfinite(est.step_between) else "not measured"
        sa = f"{np.degrees(est.step_along):.5g}" if np.isfinite(est.step_along) else "not measured"
        L.append(f"                 spacing of the returns: {sb} deg between lines, {sa} deg along them "
                 f"(declared step {np.degrees(scan.step_az):.5g} x {np.degrees(scan.step_zen):.5g})")
    elif est is not None:
        L.append("  scan lines     not analysed: " + est.reason)
    for w in info["warnings"]:
        L.append("  " + w.split(": ", 1)[-1] + "   <-- check")
    lo, hi = scan.xyz.min(axis=0), scan.xyz.max(axis=0)
    L.append(f"  extent         x {lo[0]:.2f}..{hi[0]:.2f}  y {lo[1]:.2f}..{hi[1]:.2f}  z {lo[2]:.2f}..{hi[2]:.2f}")
    L.append(f"  range          {rng.min():.2f}..{rng.max():.2f} m (median {np.median(rng):.2f})")
    zd = np.degrees(zen)
    L.append(f"  zenith angles  {zd.min():.1f}..{zd.max():.1f} deg seen; declared field of view "
             f"{np.degrees(scan.zenith_min):.1f}..{np.degrees(scan.zenith_max):.1f}")
    mid = np.arctan2(np.mean(np.sin(az)), np.mean(np.cos(az)))
    a = np.degrees(np.mod(az - mid + np.pi, 2 * np.pi) - np.pi)
    L.append(f"  bearings       {np.mod(np.degrees(mid) + a.min(), 360):.1f}..{np.mod(np.degrees(mid) + a.max(), 360):.1f} deg "
             f"clockwise from +y seen (span {a.max() - a.min():.1f} deg)")
    if scan.left_edge is None:
        L.append("  declared sector full circle")
    else:
        o = scan.position
        b_l = np.degrees(np.arctan2(scan.left_edge[0] - o[0], scan.left_edge[1] - o[1]))
        b_r = np.degrees(np.arctan2(scan.right_edge[0] - o[0], scan.right_edge[1] - o[1]))
        width = np.mod(b_r - b_l, 360.0)
        frac = float(np.mean(geo.in_sector_about(scan.xyz - o, o, scan.left_edge, scan.right_edge, info["axis"])))
        L.append(f"  declared sector {np.mod(b_l, 360):.1f}..{np.mod(b_r, 360):.1f} deg ({width if width > 0 else 360.0:.1f} deg wide); "
                 f"{100 * frac:.1f}% of returns inside" + ("   <-- check left_edge and right_edge" if frac < 0.99 else ""))
    outz = float(np.mean((zen < scan.zenith_min) | (zen > scan.zenith_max)))
    if outz > 0.001:
        L.append(f"  {100 * outz:.1f}% of returns lie outside the declared zenith range   <-- check zenith_range_deg")
    ratio, nwin = density_check(scan)
    if nwin:
        msg = f"  densest patches hold {100 * ratio:.0f}% of the pulses the declared step allows"
        if ratio > 1.15:
            msg += ("   <-- more returns than pulses: the step is finer than declared, the pattern is wrong, "
                    "or there are several returns per pulse")
        elif ratio < 0.6:
            msg += "   (no fully hit surface found, or the real step is coarser than declared)"
        L.append(msg)
    return "\n".join(L)
