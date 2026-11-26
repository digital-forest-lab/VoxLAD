"""Counting engine that reproduces the MATLAB algorithm test for test.

For each voxel and scan the MATLAB script loops over every return and applies,
in this order (voxlad_v7.m lines 658-734; check_behind.m):

  1. skip the return if it lies beyond the voxel's farthest vertex along x, y
     or z;
  2. skip it if it is farther than one voxel side from the scanner-centre line;
  3. skip it if it is on the wrong side of any of the six planes through the
     scanner and two consecutive vertices of the voxel's outline;
  4. count it as a hit (leaf, wood or noise) if it is strictly inside the
     voxel, otherwise as a blocked pulse.

The arithmetic below keeps the operation order of the MATLAB statements in
double precision, so integer counts can be compared one to one with a MATLAB
or Octave run.  It is only meant for that purpose: its cost grows with
(number of voxels) x (number of returns).
"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

# outline (six vertices, closed) for the eight cases of the MATLAB code,
# 0-based vertex numbers; rows are cases 1..8 in the MATLAB order
_OUTLINE = np.array([[3, 2, 1, 5, 4, 7, 3],
                     [0, 1, 5, 6, 7, 3, 0],
                     [7, 6, 2, 1, 0, 4, 7],
                     [3, 2, 6, 5, 4, 0, 3],
                     [0, 1, 2, 6, 7, 4, 0],
                     [7, 6, 5, 1, 0, 3, 7],
                     [4, 5, 6, 2, 3, 0, 4],
                     [4, 5, 1, 2, 3, 7, 4]], dtype=np.int64)
_FARTHEST = np.array([0, 4, 3, 1, 5, 2, 7, 6], dtype=np.int64)
_CLOSEST = np.array([6, 2, 5, 7, 3, 4, 1, 0], dtype=np.int64)


def _impl(px, py, pz, ptype, lo, vs, tls, outline, farthest, closest, out):
    x_tls = tls[0]
    y_tls = tls[1]
    z_tls = tls[2]
    npts = px.shape[0]
    for v in prange(lo.shape[0]):
        x_min = lo[v, 0]
        y_min = lo[v, 1]
        z_min = lo[v, 2]
        x_max = x_min + vs
        y_max = y_min + vs
        z_max = z_min + vs
        cx_ = x_min + (vs / 2)
        cy_ = y_min + (vs / 2)
        cz_ = z_min + (vs / 2)
        vert = np.empty((8, 3))
        vert[0, 0] = x_min; vert[0, 1] = y_min; vert[0, 2] = z_min
        vert[1, 0] = x_min; vert[1, 1] = y_min; vert[1, 2] = z_max
        vert[2, 0] = x_max; vert[2, 1] = y_min; vert[2, 2] = z_max
        vert[3, 0] = x_max; vert[3, 1] = y_min; vert[3, 2] = z_min
        vert[4, 0] = x_min; vert[4, 1] = y_max; vert[4, 2] = z_min
        vert[5, 0] = x_min; vert[5, 1] = y_max; vert[5, 2] = z_max
        vert[6, 0] = x_max; vert[6, 1] = y_max; vert[6, 2] = z_max
        vert[7, 0] = x_max; vert[7, 1] = y_max; vert[7, 2] = z_min
        # which of each pair of faces is nearer (MATLAB variables cx, cy, cz)
        cx = -1 if abs(y_tls - vert[0, 1]) < abs(y_tls - vert[4, 1]) else 1
        cy = -1 if abs(x_tls - vert[0, 0]) < abs(x_tls - vert[3, 0]) else 1
        cz = -1 if abs(z_tls - vert[0, 2]) < abs(z_tls - vert[1, 2]) else 1
        if cx == 1 and cy == 1 and cz == 1:
            cas = 0
        elif cx == -1 and cy == 1 and cz == 1:
            cas = 1
        elif cx == 1 and cy == -1 and cz == 1:
            cas = 2
        elif cx == 1 and cy == 1 and cz == -1:
            cas = 3
        elif cx == -1 and cy == 1 and cz == -1:
            cas = 4
        elif cx == 1 and cy == -1 and cz == -1:
            cas = 5
        elif cx == -1 and cy == -1 and cz == 1:
            cas = 6
        else:
            cas = 7
        fp = farthest[cas]
        cp = closest[cas]
        far_x = vert[fp, 0]
        far_y = vert[fp, 1]
        far_z = vert[fp, 2]
        delta_x = vert[fp, 0] - vert[cp, 0]
        delta_y = vert[fp, 1] - vert[cp, 1]
        delta_z = vert[fp, 2] - vert[cp, 2]
        pos_x = (abs(delta_x) - delta_x) == 0
        pos_y = (abs(delta_y) - delta_y) == 0
        pos_z = (abs(delta_z) - delta_z) == 0
        # the six planes (they do not depend on the return)
        A = np.empty(6)
        B = np.empty(6)
        C = np.empty(6)
        D = np.empty(6)
        nrm = np.empty(6)
        dc = np.empty(6)
        for i in range(6):
            a1 = vert[outline[cas, i], 0] - x_tls
            a2 = vert[outline[cas, i + 1], 0] - x_tls
            b1 = vert[outline[cas, i], 1] - y_tls
            b2 = vert[outline[cas, i + 1], 1] - y_tls
            c1 = vert[outline[cas, i], 2] - z_tls
            c2 = vert[outline[cas, i + 1], 2] - z_tls
            A[i] = b1 * c2 - b2 * c1
            B[i] = -(a1 * c2 - a2 * c1)
            C[i] = a1 * b2 - a2 * b1
            D[i] = -A[i] * x_tls - B[i] * y_tls - C[i] * z_tls
            nrm[i] = np.sqrt(A[i] ** 2 + B[i] ** 2 + C[i] ** 2)
            dc[i] = (A[i] * cx_ + B[i] * cy_ + C[i] * cz_ + D[i]) / nrm[i]
        vx = cx_ - x_tls
        vy = cy_ - y_tls
        vz = cz_ - z_tls
        vnorm = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        nb = 0
        nl = 0
        nw = 0
        nn = 0
        for p in range(npts):
            x = px[p]
            y = py[p]
            z = pz[p]
            if pos_x:
                if x > far_x:
                    continue
            elif x < far_x:
                continue
            if pos_y:
                if y > far_y:
                    continue
            elif y < far_y:
                continue
            if pos_z:
                if z > far_z:
                    continue
            elif z < far_z:
                continue
            ux = x - x_tls
            uy = y - y_tls
            uz = z - z_tls
            d_p = np.sqrt((uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2) / vnorm
            if d_p > vs:
                continue
            marker = False
            for i in range(6):
                dist_pt = (A[i] * x + B[i] * y + C[i] * z + D[i]) / nrm[i]
                if (dc[i] > 0 and dist_pt < 0) or (dc[i] < 0 and dist_pt > 0):
                    marker = True
                    break
            if marker:
                continue
            if x > x_min and x < x_max and y > y_min and y < y_max and z > z_min and z < z_max:
                t = ptype[p]
                if t < 2:
                    nl += 1
                if t == 2:
                    nw += 1
                if t > 2:
                    nn += 1
            else:
                nb += 1
        out[v, 0] = nb
        out[v, 1] = nl
        out[v, 2] = nw
        out[v, 3] = nn


if HAVE_NUMBA:
    _impl_jit = njit(parallel=True, cache=True, nogil=True)(_impl)


def count_matlab(xyz, cls, lo, voxel_size, scanner):
    """Counts (N_B, N_I leaf, N_I wood, N_I noise) per voxel, MATLAB rules.

    xyz : (n, 3) returns of one scan;  cls : class codes 1 leaf, 2 wood, 3 noise
    lo  : (m, 3) lower corners of the voxels to evaluate
    Returns an (m, 4) int64 array.
    """
    if not HAVE_NUMBA:
        raise RuntimeError(
            "the MATLAB-compatible engine needs numba (pip install numba); "
            "the default engine does not")
    xyz = np.asarray(xyz, dtype=np.float64)
    out = np.zeros((np.shape(lo)[0], 4), dtype=np.int64)
    _impl_jit(np.ascontiguousarray(xyz[:, 0]), np.ascontiguousarray(xyz[:, 1]),
              np.ascontiguousarray(xyz[:, 2]), np.ascontiguousarray(cls, dtype=np.int64),
              np.ascontiguousarray(lo, dtype=np.float64), float(voxel_size),
              np.asarray(scanner, dtype=np.float64), _OUTLINE, _FARTHEST, _CLOSEST, out)
    return out
