"""voxlad_v7.m stores coordinates in single precision.  Compare its Octave run
on the 1:400 thinned scans with the same run using double precision storage."""
import os

import numpy as np

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "octave", "reference")


def load(tag):
    d = np.loadtxt(os.path.join(REF, tag + "_dbg.txt"))
    return {(round(r[0], 2), round(r[1], 2), round(r[2], 2), int(r[3])): r for r in d}


s, d = load("thin400_single"), load("thin400")
common = sorted(set(s) & set(d))
vox = lambda t: {k[:3] for k in t}
print(f"voxels evaluated: single {len(vox(s))}, double {len(vox(d))}; "
      f"only in single {len(vox(s) - vox(d))}, only in double {len(vox(d) - vox(s))}")
ys, yd = sorted({k[1] for k in s}), sorted({k[1] for k in d})
print(f"y of the last voxel layer: single {ys[-1]}, double {yd[-1]} "
      f"({len({k[:3] for k in d if k[1] == yd[-1]})} voxels in that layer)")
for name, c in (("N_T", 8), ("N_B", 9), ("N_I leaf", 10), ("N_I wood", 11), ("N_I noise", 12)):
    a = np.array([s[k][c] for k in common])
    b = np.array([d[k][c] for k in common])
    print(f"{name:9s} differs in {int((a != b).sum())} of {len(common)} voxel-scans, by at most {np.abs(a - b).max():.0f}")
for name, c in (("A_proj", 4), ("G", 5)):
    a = np.array([s[k][c] for k in common])
    b = np.array([d[k][c] for k in common])
    print(f"{name:9s} largest relative difference {np.max(np.abs(a - b) / np.abs(b)):.1e}")
