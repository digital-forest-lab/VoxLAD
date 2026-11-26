"""Fast counting engine: follow each pulse through the voxel grid.

For every return P of a scan taken from position O, the pulse is known to
have travelled from O to P.  Every voxel that the straight line from O through
P crosses *beyond* P is a voxel the pulse was heading for but never reached:
P is one of that voxel's "blocked" pulses (N_B in the paper).  So instead of
testing every return against every voxel (the MATLAB algorithm), each return
is walked once through the grid with a 3-D digital differential analyser
(Amanatides & Woo 1987), adding one to N_B of each voxel on the way.  The
cost is proportional to (number of returns) x (voxels along a ray), and does
not depend on how many voxels are being evaluated.
"""
from __future__ import annotations

import numpy as np

try:  # numba is optional; a NumPy fallback is used without it
    import numba
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False

_BIG = 1.0e300


def point_voxel_index(xyz, origin, h, shape):
    """Flat voxel index of each point (-1 outside the grid).

    A point on a face shared by two voxels belongs to the voxel on the upper
    side: voxel i spans [origin + i*h, origin + (i+1)*h).
    """
    ijk = np.floor((np.asarray(xyz, dtype=float) - np.asarray(origin, dtype=float)) / h).astype(np.int64)
    nx, ny, nz = shape
    inside = ((ijk[:, 0] >= 0) & (ijk[:, 0] < nx) & (ijk[:, 1] >= 0) & (ijk[:, 1] < ny)
              & (ijk[:, 2] >= 0) & (ijk[:, 2] < nz))
    flat = (ijk[:, 0] * ny + ijk[:, 1]) * nz + ijk[:, 2]
    return np.where(inside, flat, -1)


if HAVE_NUMBA:

    @njit(cache=True, nogil=True)
    def _walk(px, py, pz, ox, oy, oz, gx, gy, gz, h, nx, ny, nz, lookup, acc):
        dx = px - ox
        dy = py - oy
        dz = pz - oz
        if dx == 0.0 and dy == 0.0 and dz == 0.0:
            return
        # voxel holding the return, if it is inside the grid
        fi = np.floor((px - gx) / h)
        fj = np.floor((py - gy) / h)
        fk = np.floor((pz - gz) / h)
        inside = (fi >= 0 and fi < nx and fj >= 0 and fj < ny and fk >= 0 and fk < nz)
        if inside:
            i = int(fi)
            j = int(fj)
            k = int(fk)
        else:
            # where does the line beyond the return enter the grid box?
            t0 = 0.0
            t1 = _BIG
            lo = gx
            hi = gx + nx * h
            if dx != 0.0:
                ta = (lo - px) / dx
                tb = (hi - px) / dx
                if ta > tb:
                    ta, tb = tb, ta
                if ta > t0:
                    t0 = ta
                if tb < t1:
                    t1 = tb
            elif px < lo or px >= hi:
                return
            lo = gy
            hi = gy + ny * h
            if dy != 0.0:
                ta = (lo - py) / dy
                tb = (hi - py) / dy
                if ta > tb:
                    ta, tb = tb, ta
                if ta > t0:
                    t0 = ta
                if tb < t1:
                    t1 = tb
            elif py < lo or py >= hi:
                return
            lo = gz
            hi = gz + nz * h
            if dz != 0.0:
                ta = (lo - pz) / dz
                tb = (hi - pz) / dz
                if ta > tb:
                    ta, tb = tb, ta
                if ta > t0:
                    t0 = ta
                if tb < t1:
                    t1 = tb
            elif pz < lo or pz >= hi:
                return
            if t0 >= t1:
                return
            # step a hair inside the box to find the entry voxel
            tm = t0 + 1e-9 * (t1 - t0)
            i = int(np.floor((px + tm * dx - gx) / h))
            j = int(np.floor((py + tm * dy - gy) / h))
            k = int(np.floor((pz + tm * dz - gz) / h))
            if i < 0:
                i = 0
            elif i >= nx:
                i = nx - 1
            if j < 0:
                j = 0
            elif j >= ny:
                j = ny - 1
            if k < 0:
                k = 0
            elif k >= nz:
                k = nz - 1
            c = lookup[(i * ny + j) * nz + k]
            if c >= 0:
                acc[c] += 1
        # set up the walk; t is measured from the return along the pulse
        if dx > 0.0:
            sx = 1
            tmx = (gx + (i + 1) * h - px) / dx
            tdx = h / dx
        elif dx < 0.0:
            sx = -1
            tmx = (gx + i * h - px) / dx
            tdx = -h / dx
        else:
            sx = 0
            tmx = _BIG
            tdx = _BIG
        if dy > 0.0:
            sy = 1
            tmy = (gy + (j + 1) * h - py) / dy
            tdy = h / dy
        elif dy < 0.0:
            sy = -1
            tmy = (gy + j * h - py) / dy
            tdy = -h / dy
        else:
            sy = 0
            tmy = _BIG
            tdy = _BIG
        if dz > 0.0:
            sz = 1
            tmz = (gz + (k + 1) * h - pz) / dz
            tdz = h / dz
        elif dz < 0.0:
            sz = -1
            tmz = (gz + k * h - pz) / dz
            tdz = -h / dz
        else:
            sz = 0
            tmz = _BIG
            tdz = _BIG
        while True:
            if tmx <= tmy and tmx <= tmz:
                i += sx
                tmx += tdx
                if i < 0 or i >= nx:
                    return
            elif tmy <= tmz:
                j += sy
                tmy += tdy
                if j < 0 or j >= ny:
                    return
            else:
                k += sz
                tmz += tdz
                if k < 0 or k >= nz:
                    return
            c = lookup[(i * ny + j) * nz + k]
            if c >= 0:
                acc[c] += 1

    @njit(parallel=True, cache=True, nogil=True)
    def _blocked_numba(xyz, ox, oy, oz, gx, gy, gz, h, nx, ny, nz, lookup, nsel, nchunks):
        n = xyz.shape[0]
        acc = np.zeros((nchunks, nsel), dtype=np.int64)
        for c in prange(nchunks):
            a = acc[c]
            for p in range(c * n // nchunks, (c + 1) * n // nchunks):
                _walk(xyz[p, 0], xyz[p, 1], xyz[p, 2], ox, oy, oz, gx, gy, gz, h, nx, ny, nz, lookup, a)
        return acc.sum(axis=0)


def _blocked_numpy(xyz, o, g0, h, shape, lookup, nsel, batch=2_000_000):
    """Same walk as _walk, all pulses advanced together one voxel per pass."""
    nx, ny, nz = shape
    dims = np.array([nx, ny, nz])
    out = np.zeros(nsel, dtype=np.int64)
    box_lo = np.asarray(g0, dtype=float)
    box_hi = box_lo + dims * h
    for s in range(0, xyz.shape[0], batch):
        p = np.asarray(xyz[s:s + batch], dtype=float)
        d = p - o
        keep = np.any(d != 0, axis=1)
        p, d = p[keep], d[keep]
        fidx = np.floor((p - box_lo) / h)
        inside = np.all((fidx >= 0) & (fidx < dims), axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            ta = (box_lo - p) / d
            tb = (box_hi - p) / d
        par = d == 0
        tlo = np.where(par, -np.inf, np.minimum(ta, tb))
        thi = np.where(par, np.inf, np.maximum(ta, tb))
        par_out = np.any(par & ((p < box_lo) | (p >= box_hi)), axis=1)
        t0 = np.maximum(tlo.max(axis=1), 0.0)
        t1 = thi.min(axis=1)
        enters = ~inside & ~par_out & (t0 < t1)
        tm = t0 + 1e-9 * (np.where(np.isfinite(t1), t1, t0) - t0)
        eidx = np.clip(np.floor((p + tm[:, None] * d - box_lo) / h), 0, dims - 1)
        idx = np.where(inside[:, None], fidx, eidx).astype(np.int64)
        act = inside | enters
        # entry voxel of pulses that start outside the grid is already "beyond"
        e = idx[enters]
        cc = lookup[(e[:, 0] * ny + e[:, 1]) * nz + e[:, 2]]
        out += np.bincount(cc[cc >= 0], minlength=nsel)
        p, d, idx = p[act], d[act], idx[act]
        step = np.sign(d).astype(np.int64)
        with np.errstate(divide="ignore", invalid="ignore"):
            nxt = box_lo + (idx + (step > 0)) * h
            tmax = np.where(step != 0, (nxt - p) / d, _BIG)
            tdel = np.where(step != 0, h / np.abs(d), _BIG)
        while idx.shape[0]:
            # same tie rule as the compiled walk: x before y before z
            ax = np.where((tmax[:, 0] <= tmax[:, 1]) & (tmax[:, 0] <= tmax[:, 2]), 0,
                          np.where(tmax[:, 1] <= tmax[:, 2], 1, 2))
            r = np.arange(idx.shape[0])
            idx[r, ax] += step[r, ax]
            tmax[r, ax] += tdel[r, ax]
            ok = (idx[r, ax] >= 0) & (idx[r, ax] < dims[ax])
            idx, step, tmax, tdel = idx[ok], step[ok], tmax[ok], tdel[ok]
            cc = lookup[(idx[:, 0] * ny + idx[:, 1]) * nz + idx[:, 2]]
            out += np.bincount(cc[cc >= 0], minlength=nsel)
    return out


def count_blocked(xyz, scanner, grid_origin, h, shape, lookup, nsel, threads=None, use_numba=None):
    """N_B for the selected voxels: number of returns lying on a line of sight
    from the scanner to the voxel, in front of it.

    lookup : int32 array of size nx*ny*nz; position of each voxel in the
             output, or -1 for voxels that are not wanted.
    """
    if use_numba is None:
        use_numba = HAVE_NUMBA
    xyz = np.ascontiguousarray(xyz, dtype=np.float64)
    o = np.asarray(scanner, dtype=float)
    g0 = np.asarray(grid_origin, dtype=float)
    nx, ny, nz = (int(v) for v in shape)
    lookup = np.ascontiguousarray(lookup, dtype=np.int32)
    if use_numba and HAVE_NUMBA:
        if threads:
            numba.set_num_threads(max(1, min(int(threads), numba.config.NUMBA_NUM_THREADS)))
        nchunks = max(1, min(numba.get_num_threads(), xyz.shape[0] // 10000 + 1))
        # keep the per-thread accumulators under about 1 GB in total
        nchunks = max(1, min(nchunks, int(1e9 // max(8 * nsel, 1)) or 1))
        return _blocked_numba(xyz, o[0], o[1], o[2], g0[0], g0[1], g0[2], float(h),
                              nx, ny, nz, lookup, int(nsel), int(nchunks))
    return _blocked_numpy(xyz, o, g0, float(h), (nx, ny, nz), lookup, int(nsel))


def count_hits(xyz, cls, grid_origin, h, shape, lookup, nsel):
    """Returns inside each selected voxel, by class (1 leaf, 2 wood, 3 noise).

    Gives an (nsel, 3) integer array.
    """
    flat = point_voxel_index(xyz, grid_origin, h, shape)
    ok = flat >= 0
    comp = np.full(flat.shape, -1, dtype=np.int64)
    comp[ok] = lookup[flat[ok]]
    out = np.zeros((nsel, 3), dtype=np.int64)
    for c in (1, 2, 3):
        m = (comp >= 0) & (cls == c)
        out[:, c - 1] = np.bincount(comp[m], minlength=nsel)
    return out
