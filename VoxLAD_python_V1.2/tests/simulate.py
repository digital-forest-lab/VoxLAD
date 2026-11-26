"""A small lidar simulator used by the tests: a box of turbid medium (leaves as
a 'gas' of small scatterers) scanned from one or more positions.

Each pulse travels in a straight line; inside the box it is intercepted with
probability 1 - exp(-LAD * G(zenith) * path), the return being placed at an
exponentially distributed depth.  Pulses that are not intercepted give no
return, as when they end in the sky.
"""
import numpy as np

from voxlad import leafangle


def tilt_matrix(tilt_deg, bearing_deg, heading_deg=0.0):
    """Rotation (scanner frame -> common frame) of a scanner whose head axis leans by ``tilt_deg``
    from the vertical towards the compass-like bearing ``bearing_deg`` (clockwise from +y)."""
    def rot(axis, ang):
        axis = np.asarray(axis, float) / np.linalg.norm(axis)
        K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
    b = np.radians(bearing_deg)
    lean = np.array([np.sin(b), np.cos(b), 0.0])                 # horizontal direction of the lean
    return rot(np.cross([0, 0, 1.0], lean), np.radians(tilt_deg)) @ rot([0, 0, 1.0], np.radians(heading_deg))


def scan_room(position, rotation, step, rng, room_lo=(-12.0, -10.0, -1.6), room_hi=(11.0, 13.0, 9.0),
              polar_range=(30.0, 130.0), head_range=(0.0, 360.0), cone_deg=0.0, origin_offset=0.0, dropout=0.2):
    """A panoramic scanner in a closed room: every pulse returns from a wall, the floor or the ceiling,
    except a random share.  The returns come line by line, as in a file kept in the order of acquisition.

    rotation       scanner frame -> common frame (the head axis is the scanner's z)
    cone_deg       half-angle by which the scan lines are cones and not planes through the axis
    origin_offset  lateral offset of the beam origin (m), along the axis of each line cone
    Returns (xyz, head axis in the common frame).
    """
    o = np.asarray(position, float)
    R = np.eye(3) if rotation is None else np.asarray(rotation, float)
    tau = np.radians(np.arange(polar_range[0], polar_range[1], np.degrees(step)))
    psi = np.radians(np.arange(head_range[0], head_range[1], np.degrees(step)))
    P, T = np.meshgrid(psi, tau, indexing="ij")                 # line by line
    P, T = P.ravel(), T.ravel()
    m = np.column_stack([np.cos(P), np.sin(P), np.zeros_like(P)])
    side = np.column_stack([-np.sin(P), np.cos(P), np.zeros_like(P)])       # z x m
    d = np.sin(np.radians(cone_deg)) * m + np.cos(np.radians(cone_deg)) * (
        np.cos(T)[:, None] * np.array([0, 0, 1.0]) + np.sin(T)[:, None] * side)
    d = d @ R.T
    org = o + origin_offset * (m @ R.T)
    lo, hi = np.asarray(room_lo, float), np.asarray(room_hi, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(d > 0, (hi - org) / d, np.where(d < 0, (lo - org) / d, np.inf)).min(axis=1)
    keep = rng.random(len(t)) >= dropout
    return (org + d * t[:, None])[keep], R[:, 2].copy()


def scan_box(position, box_lo, box_hi, lad, g_table, step_az, step_zen, rng, pattern="spherical", rotation=None):
    """Returns (xyz of returns, number of pulses aimed at the box).

    ``rotation`` (scanner frame -> common frame) tilts the scanner: the pulses are then equally
    spaced in azimuth and polar angle about its own axis, while the box and the leaf angles keep
    to the common frame.  The returns come line by line.
    """
    o = np.asarray(position, float)
    lo, hi = np.asarray(box_lo, float), np.asarray(box_hi, float)
    R = np.eye(3) if rotation is None else np.asarray(rotation, float)
    # angular extent of the box: sample its surface (the extreme zenith angles are
    # reached on edges and faces, not only at corners)
    t = np.linspace(0.0, 1.0, 41)
    u, v = np.meshgrid(t, t, indexing="ij")
    faces = []
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        for side in (lo[a], hi[a]):
            p = np.empty(u.shape + (3,))
            p[..., a] = side
            p[..., b] = lo[b] + u * (hi[b] - lo[b])
            p[..., c] = lo[c] + v * (hi[c] - lo[c])
            faces.append(p.reshape(-1, 3))
    surf = (np.vstack(faces) - o) @ R                           # the box seen in the scanner's own frame
    az = np.arctan2(surf[:, 0], surf[:, 1])
    zen = np.arccos(surf[:, 2] / np.linalg.norm(surf, axis=1))
    pad = 8
    a = np.arange(np.floor(az.min() / step_az) - pad, np.ceil(az.max() / step_az) + pad + 1) * step_az
    z = np.arange(np.floor(zen.min() / step_zen) - pad, np.ceil(zen.max() / step_zen) + pad + 1) * step_zen
    A, Z = np.meshgrid(a, z, indexing="ij")
    A, Z = A.ravel(), Z.ravel()
    if pattern == "uniform":
        # same spacing in both directions: widen the azimuth step by 1/sin(zenith) row by row
        rows = []
        for zz in z:
            st = step_az / np.sin(zz)
            aa = np.arange(np.floor(az.min() / st) - pad, np.ceil(az.max() / st) + pad + 1) * st
            rows.append(np.column_stack([aa, np.full(aa.size, zz)]))
        AZ = np.vstack(rows)
        A, Z = AZ[:, 0], AZ[:, 1]
    d = np.column_stack([np.sin(Z) * np.sin(A), np.sin(Z) * np.cos(A), np.cos(Z)]) @ R.T
    Z = np.arccos(np.clip(d[:, 2], -1.0, 1.0))                  # the leaves see the true zenith angle
    with np.errstate(divide="ignore", invalid="ignore"):
        ta = (lo - o) / d
        tb = (hi - o) / d
    t0 = np.minimum(ta, tb).max(axis=1)
    t1 = np.maximum(ta, tb).min(axis=1)
    hit = (t1 > t0) & (t1 > 0)
    t0 = np.maximum(t0[hit], 0.0)
    path = t1[hit] - t0
    k = lad * leafangle.g_function(Z[hit], g_table)
    s = rng.exponential(1.0, size=path.size) / k
    got = s < path
    xyz = o + d[hit][got] * (t0[got] + s[got])[:, None]
    return xyz, int(hit.sum())
