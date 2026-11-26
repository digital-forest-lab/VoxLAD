"""Timings on the two example scans (1.8 million returns, 0.3 m voxels)."""
import os
import time

import numpy as np

from _common import EXAMPLE, example_scans, voxlad
from voxlad import traverse
from voxlad.model import Settings, compute

t = time.perf_counter()
scans, settings = example_scans()
t_read = time.perf_counter() - t
n = sum(s.xyz.shape[0] for s in scans)
print(f"reading the two .asc files ({n} returns): {t_read:.2f} s")
compute([voxlad.Scan(s.xyz[::500], s.cls[::500], s.position, s.pattern, s.step_az) for s in scans],
        Settings(voxel_size=0.3, z_range=(-1, 12)))                    # compile
try:
    import numba
    threads = numba.get_num_threads()
except ImportError:
    threads = 1
print(f"threads: {threads}; cores: {os.cpu_count()}")


def best(fn, rep=5):
    ts = []
    for _ in range(rep):
        t = time.perf_counter()
        r = fn()
        ts.append(time.perf_counter() - t)
    return min(ts), r


for h in (0.3, 0.1):
    s = Settings(voxel_size=h, leaf_angle="erectophile", z_range=(-1, 12))
    dt, r = best(lambda: compute(scans, s))
    print(f"default engine, {h} m voxels: {dt:.2f} s for {r.lad.size} voxels "
          f"(counting {r.timings['counting']:.2f} s)")
s = Settings(voxel_size=0.3, leaf_angle="erectophile", z_range=(-1, 12))
r = compute(scans, s)
lookup = np.full(int(np.prod(r.shape)), -1, np.int32)
flat = (r.ijk[:, 0] * r.shape[1] + r.ijk[:, 1]) * r.shape[2] + r.ijk[:, 2]
lookup[flat] = np.arange(flat.size)
dt, _ = best(lambda: [traverse.count_blocked(sc.xyz, sc.position, r.origin, 0.3, r.shape, lookup, flat.size,
                                             use_numba=False) for sc in scans], rep=2)
print(f"counting without Numba (NumPy fallback), 0.3 m: {dt:.2f} s")
if traverse.HAVE_NUMBA:
    t = time.perf_counter()
    rm = compute(scans, Settings(voxel_size=0.3, leaf_angle="erectophile", z_range=(-1, 12), compat="matlab"))
    dt = time.perf_counter() - t
    nused = sum(int(((sc.xyz[:, 2] < 12) & (sc.xyz[:, 2] > -1)).sum()) for sc in scans)   # after the height filter
    print(f"MATLAB algorithm compiled with Numba (compat = matlab): {dt:.1f} s for {rm.lad.size} voxels, "
          f"i.e. {rm.lad.size} x {nused} = {rm.lad.size * nused:.2e} return-voxel tests")
# a larger job: the scans repeated with small shifts
big = []
rng = np.random.default_rng(0)
for sc in scans:
    xyz = np.vstack([sc.xyz + rng.normal(0, 0.002, sc.xyz.shape) for _ in range(10)])
    big.append(voxlad.Scan(xyz, np.tile(sc.cls, 10), sc.position, sc.pattern, sc.step_az / np.sqrt(10)))
dt, r = best(lambda: compute(big, Settings(voxel_size=0.1, leaf_angle="erectophile", z_range=(-1, 12))), rep=2)
print(f"ten times the returns ({10 * n}), 0.1 m voxels: {dt:.2f} s for {r.lad.size} voxels")
