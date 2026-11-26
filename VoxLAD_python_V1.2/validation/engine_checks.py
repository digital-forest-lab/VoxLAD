"""Check the fast (traversal) engine on the two example scans.

1. against an independent brute-force test (slab-method ray/box intersection
   in NumPy) on 400 random voxel-scans;
2. NumPy fallback against the compiled walk;
3. against the MATLAB-rule engine, separating the voxels for which the scanner
   lies inside the voxel's own x, y or z extent (where the MATLAB outline of
   the voxel is wrong) from the others.

The grid is shifted by 0.05 mm so that no millimetre-quantised return lies
exactly on a voxel face; the effect of face returns is reported separately on
the unshifted grid.
"""
import numpy as np

from _common import example_scans, voxlad
from voxlad import traverse
from voxlad.model import Settings, compute

scans, _ = example_scans()
h, e = 0.3, 5e-5
base = dict(voxel_size=h, leaf_angle="erectophile")
shifted = [-5 + e, 17 + e, -1 + e, 4 + e, 28 + e, 7 + e]
rt = compute(scans, Settings(bounds=shifted, **base))
rm = compute(scans, Settings(bounds=shifted, engine="matlab", **base))
print(f"shifted grid: {rt.lad.size} voxels, {rt.per_scan['nb'].size} voxel-scans")


def brute(P, O, l):
    d = P - O
    with np.errstate(divide="ignore", invalid="ignore"):
        ta, tb = (l - O) / d, (l + h - O) / d
    t0 = np.minimum(ta, tb).max(axis=1)             # the return is at t = 1 on its own line of sight
    t1 = np.maximum(ta, tb).min(axis=1)
    inside = np.all((P >= l) & (P < l + h), axis=1)
    return int(np.sum((t0 < t1) & (t0 >= 1) & ~inside)), int(inside.sum())


rng = np.random.default_rng(3)
bad_b = bad_i = 0
for _ in range(400):
    v, s = int(rng.integers(rt.lad.size)), int(rng.integers(2))
    nb, ni = brute(scans[s].xyz, scans[s].position, rt.lo[v])
    bad_b += nb != rt.per_scan["nb"][v, s]
    bad_i += ni != rt.per_scan["ni_leaf"][v, s] + rt.per_scan["ni_wood"][v, s] + rt.per_scan["ni_noise"][v, s]
print(f"1. brute force vs traversal, 400 random voxel-scans: N_B differs in {bad_b}, N_I differs in {bad_i}")

lookup = np.full(int(np.prod(rt.shape)), -1, np.int32)
flat = (rt.ijk[:, 0] * rt.shape[1] + rt.ijk[:, 1]) * rt.shape[2] + rt.ijk[:, 2]
lookup[flat] = np.arange(flat.size)
same = all(np.array_equal(traverse.count_blocked(sc.xyz, sc.position, rt.origin, h, rt.shape, lookup, flat.size,
                                                 use_numba=False), rt.per_scan["nb"][:, j])
           for j, sc in enumerate(scans))
print(f"2. NumPy fallback identical to the compiled walk for both scans: {same}")

d = rm.per_scan["nb"] - rt.per_scan["nb"]
inslab = np.stack([np.any((sc.position >= rt.lo) & (sc.position <= rt.lo + h), axis=1) for sc in scans], axis=1)
print(f"3. MATLAB-rule engine minus traversal; scanner inside the voxel's x, y or z extent for "
      f"{inslab.sum()} of {inslab.size} voxel-scans")
for name, m in (("scanner outside", ~inslab), ("scanner inside ", inslab)):
    tot = rt.per_scan["nb"][m].sum()
    print(f"   {name}: N_B identical in {100 * np.mean(d[m] == 0):.2f}% of voxel-scans, largest difference "
          f"{np.abs(d[m]).max()}, sum {d[m].sum():+d} of {tot} ({100 * d[m].sum() / tot:+.3f}%)")
rel = d[inslab] / np.maximum(rt.per_scan["nb"][inslab], 1)
print(f"   scanner inside: median relative difference {100 * np.median(rel):+.2f}%")
w = np.unravel_index(np.argmin(d), d.shape)
print(f"   largest: voxel at {np.round(rt.lo[w[0]], 2).tolist()}, scan {w[1] + 1}: "
      f"{rm.per_scan['nb'][w]} (MATLAB rules) vs {rt.per_scan['nb'][w]}")

# returns lying exactly on voxel faces: unshifted grid
nominal = [-5, 17, -1, 4, 28, 7]
a = compute(scans, Settings(bounds=nominal, **base))
b = compute(scans, Settings(bounds=nominal, engine="matlab", **base))
print(f"4. unshifted grid, leaf returns counted inside voxels: MATLAB rules {b.per_scan['ni_leaf'].sum()}, "
      f"half-open voxels {a.per_scan['ni_leaf'].sum()} "
      f"({100 * (b.per_scan['ni_leaf'].sum() / a.per_scan['ni_leaf'].sum() - 1):+.2f}%)")
