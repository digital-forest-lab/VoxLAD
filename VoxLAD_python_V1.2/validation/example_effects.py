"""Effect of the corrections and of the voxel size on the example tree.

1. MATLAB-compatible mode against the default, on the same grid, voxel by voxel.
2. Total leaf area for several voxel sizes and path-length treatments.
"""
import numpy as np

from _common import example_scans, voxlad  # noqa: F401
from voxlad.model import Settings, compute

scans, _ = example_scans()
base = dict(leaf_angle="erectophile")
m = compute(scans, Settings(voxel_size=0.3, z_range=(-1, 12), compat="matlab", **base))
d = compute(scans, Settings(voxel_size=0.3, bounds=[-5, 17, -1, 4, 28, 7], **base))     # the grid of the MATLAB run
tm = np.nansum(m.lad) * 0.3 ** 3
print(f"1. MATLAB mode: {m.lad.size} voxels, {tm:.3f} m2; default on the same grid: {d.lad.size} voxels, "
      f"{d.total_leaf_area():.3f} m2 ({100 * (d.total_leaf_area() / tm - 1):+.2f}%)")
km = {tuple(k): i for i, k in enumerate(m.ijk.tolist())}
pairs = [(km[tuple(k)], i) for i, k in enumerate(d.ijk.tolist()) if tuple(k) in km]
pm, pd = np.array([a for a, _ in pairs]), np.array([b for _, b in pairs])
a, b = m.lad[pm], d.lad[pd]
ok = np.isfinite(a) & np.isfinite(b)
dd = (b - a)[ok]
print(f"   common voxels {len(pm)}; same scan retained in {100 * np.mean(m.best_scan[pm] == d.best_scan[pd]):.2f}%")
print(f"   LAD default - MATLAB: mean {dd.mean():+.4f}, rms {np.sqrt((dd ** 2).mean()):.4f} m2/m3 (mean LAD {a[ok].mean():.3f})")
for thr in (0.01, 0.05, 0.10):
    big = (np.abs(dd) > thr * np.maximum(a[ok], 1e-9)) & (np.abs(dd) > 0.005)
    print(f"   voxels differing by more than {100 * thr:.0f}%: {big.sum()} of {ok.sum()} ({100 * big.mean():.1f}%)")

print("2. total leaf area (m2) by voxel size and path-length treatment; unresolved voxels in brackets")
print(f"{'h (m)':>6s} {'voxels':>7s} {'paper':>10s} {'matlab':>12s} {'exact':>8s} {'none':>8s}")
for h in (0.1, 0.15, 0.2, 0.3, 0.5, 1.0):
    row = []
    for fit in ("paper", "matlab", "exact", "none"):
        r = compute(scans, Settings(voxel_size=h, z_range=(-1, 12), jensen_fit=fit, **base))
        row.append((r.total_leaf_area(), int((r.flag == voxlad.FLAG_UNRESOLVED).sum()), r.lad.size))
    print(f"{h:6.2f} {row[0][2]:7d} {row[0][0]:7.2f} ({row[0][1]}) {row[1][0]:7.2f} ({row[1][1]:2d}) {row[2][0]:8.2f} {row[3][0]:8.2f}")
