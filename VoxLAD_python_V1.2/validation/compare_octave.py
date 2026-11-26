"""Compare the Python MATLAB-compatible computations with the original
voxlad_v7.m run in GNU Octave.

    python compare_octave.py thin400      # scans thinned 1:400, whole grid
    python compare_octave.py thin100      # scans thinned 1:100, whole grid
    python compare_octave.py full         # full scans, 24 voxels
    python compare_octave.py all

The Octave outputs are in octave/reference (octave/make_reference.py rebuilds
the scripts that produced them).  For every voxel and scan of the Octave run,
the Python functions are called with the same voxel corner and the integer
counts, the geometry and the printed results are compared.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from voxlad import geometry as geo, lad as ladm, leafangle as la, legacy   # noqa: E402

REF = os.path.join(HERE, "octave", "reference")
SCAN_DIR = os.path.join(HERE, "..", "examples", "shea_tree")
VS = 0.3                                   # settings of voxlad_v7.m as distributed
G_TABLE = la.leaf_angle_table("erectophile")
RUNS = {"thin400": 400, "thin100": 100, "full": 1}


def read_asc(name, k):
    """Scan thinned 1:k as octave/make_reference.py writes it, after v7's height filter."""
    with open(os.path.join(SCAN_DIR, name + ".asc")) as fh:
        tls = np.array(fh.readline().split(), float)
        fh.readline(); fh.readline()
        unit, a = fh.readline().split()
        fh.readline()
        d = np.loadtxt(fh)
    alpha = float(a) * (math.pi / 180 if unit == "d" else 1.0)
    if k > 1:
        d = d[::k]
        alpha = float("%.8f" % (alpha * math.sqrt(k)))
    keep = (d[:, 2] < 12) & (d[:, 2] > -1)
    return dict(tls=tls, alpha=alpha, xyz=d[keep, :3], cls=d[keep, 3].astype(int))


def octave_axis(b, bmax, vs, split=False):
    """Loop values as Octave builds them: base + k*increment.  v7 runs x in two
    halves (its parfor), the second starting from bbox_x_min+(floor(n/2)+1)*vox_size."""
    n = int(np.floor((bmax - b) / vs + 1))
    if not split:
        return np.array([b + k * vs for k in range(n + 2)])
    half = n // 2
    lim = b + (half + 1) * vs
    return np.array([b + k * vs for k in range(half + 1)] + [lim + k * vs for k in range(n + 2)])


def fmt(v, f):
    if np.isnan(v):
        return "NaN"
    if np.isinf(v):
        return "-Inf" if v < 0 else "Inf"
    return f % v


def compare(tag):
    k = RUNS[tag]
    scans = [read_asc("a1s1_t10", k), read_asc("a1s2_t10", k)]
    allxyz = np.vstack([s["xyz"] for s in scans])
    bmin, bmax = np.floor(allxyz.min(0)), np.ceil(allxyz.max(0))
    ax = [octave_axis(bmin[0], bmax[0], VS, True), octave_axis(bmin[1], bmax[1], VS), octave_axis(bmin[2], bmax[2], VS)]
    dbg = np.loadtxt(os.path.join(REF, tag + "_dbg.txt"), ndmin=2)
    print(f"{tag}: {dbg.shape[0]} voxel-scan records of the Octave run; {sum(len(s['xyz']) for s in scans)} returns")
    lo = np.empty((dbg.shape[0], 3))
    for a in range(3):
        j = np.abs(dbg[:, a][:, None] - ax[a][None, :]).argmin(axis=1)
        assert np.all(np.abs(ax[a][j] - dbg[:, a]) < 1e-5)
        lo[:, a] = ax[a][j]
    scan_id = dbg[:, 3].astype(int)
    worst, mismatch, printed = {}, {}, {}
    for s in (1, 2):
        m = scan_id == s
        sc, L, ref = scans[s - 1], lo[m], dbg[m]
        centre = L + VS / 2
        cnt = legacy.count_matlab(sc["xyz"], sc["cls"], L, VS, sc["tls"])
        aproj, _ = geo.projected_area(L, VS, sc["tls"], face_rule="matlab")
        theta = geo.view_zenith(centre, sc["tls"])
        G = la.g_function(theta, G_TABLE, folded=False)
        dbar = VS ** 3 / aproj
        nt = geo.pulse_count(aproj, centre, sc["tls"], "spherical", sc["alpha"], sc["alpha"])   # TLS_type = 2
        for name, mine, col in (("A_proj", aproj, 4), ("G", G.real, 5), ("d_bar", dbar, 7), ("theta", theta, 13)):
            worst[name] = max(worst.get(name, 0.0), float(np.max(np.abs(mine - ref[:, col]) / np.abs(ref[:, col]))))
        worst["Im(G)"] = max(worst.get("Im(G)", 0.0), float(np.max(np.abs(G.imag - ref[:, 6]))))
        for name, mine, col in (("N_T", nt, 8), ("N_B", cnt[:, 0], 9), ("N_I leaf", cnt[:, 1], 10),
                                ("N_I wood", cnt[:, 2], 11), ("N_I noise", cnt[:, 3], 12)):
            mismatch[name] = mismatch.get(name, 0) + int(np.sum(mine != ref[:, col]))
        ng = nt - cnt.sum(axis=1)
        lad = ladm.invert_lad_matlab(dbar, ladm.jensen_delta(VS, dbar, "matlab"), G, ng, cnt[:, 1] * 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            pe = (nt - cnt[:, 0]) / nt * 100
        for i, row in enumerate(ref):
            printed[(round(row[0], 3), round(row[1], 3), round(row[2], 3), s)] = (lad[i], pe[i], nt[i] - cnt[i, 0])
    for name, v in mismatch.items():
        print(f"  {name:10s} differs in {v} records")
    for name, v in worst.items():
        print(f"  {name:10s} largest {'absolute' if name == 'Im(G)' else 'relative'} difference {v:.1e}")
    nline = nsame = nfield = nfsame = 0
    for part in ("_out_1.txt", "_out_2.txt"):
        for line in open(os.path.join(REF, tag + part)):
            t = line.split()
            if len(t) != 8:
                continue                           # the 'here' / 'there' debug lines of v7
            key = tuple(round(float(v), 3) for v in t[:3])
            if key + (1,) not in printed and key + (2,) not in printed:
                continue
            mine, ent = [], []
            for s in (1, 2):
                if key + (s,) in printed:
                    l_, p_, e_ = printed[key + (s,)]
                    mine += [fmt(l_, "%.2f"), fmt(p_, "%.0f")]
                    ent.append(e_)
                else:                               # scan skipped by v7: it prints NaN NaN
                    mine += ["NaN", "NaN"]
                    ent.append(0.0)
            mine.append(str(int(np.argmax(ent)) + 1))
            nline += 1
            nsame += mine == t[3:]
            nfield += len(mine)
            nfsame += sum(a == b for a, b in zip(mine, t[3:]))
    print(f"  result lines identical: {nsame} of {nline}; printed values identical: {nfsame} of {nfield}")
    return sum(mismatch.values()) == 0 and nsame == nline


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    tags = list(RUNS) if which == "all" else [which]
    ok = all([compare(t) for t in tags])
    print("ALL IDENTICAL" if ok else "DIFFERENCES FOUND")
    sys.exit(0 if ok else 1)
