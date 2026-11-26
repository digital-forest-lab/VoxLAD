"""Summary text and result files."""
from __future__ import annotations

import os

import numpy as np

from . import __version__ as _v  # noqa: F401  (set in __init__)
from . import io
from . import lad as ladm
from .model import FLAG_NAMES, FLAG_OK, FLAG_LOW_SAMPLING, FLAG_UNRESOLVED


def vertical_profile(res):
    """Leaf area (m2) per horizontal layer of voxels: (z_low, leaf_area, n_voxels, n_without_value)."""
    k = res.ijk[:, 2]
    la = np.where(np.isfinite(res.lad), res.lad, 0.0) * res.voxel_size ** 3
    nz = res.shape[2]
    area = np.bincount(k, weights=la, minlength=nz)
    nvox = np.bincount(k, minlength=nz)
    nnan = np.bincount(k, weights=(~np.isfinite(res.lad)).astype(float), minlength=nz).astype(int)
    zlow = res.origin[2] + np.arange(nz) * res.voxel_size
    return zlow, area, nvox, nnan


def summary_text(res):
    s, opt = res.settings, res.options
    h = res.voxel_size
    L = []
    add = L.append
    add(f"VoxLAD {_v} summary")
    add("=" * 60)
    add(f"voxel size            {h:g} m")
    add(f"grid                  {res.shape[0]} x {res.shape[1]} x {res.shape[2]} voxels, "
        f"origin ({res.origin[0]:g}, {res.origin[1]:g}, {res.origin[2]:g})")
    add(f"voxels evaluated      {res.lad.size} ({'all voxels' if s.voxels == 'all' else 'voxels holding a leaf return'})")
    la_name = s.leaf_angle if isinstance(s.leaf_angle, str) else f"{res.g_table.size} user classes"
    add(f"leaf angle            {la_name}")
    add(f"beam correction       {res.beam_factor:g}")
    add(f"settings              compat={s.compat}, engine={opt['engine']}, jensen_fit={opt['jensen_fit']}, "
        f"unresolved={opt['unresolved']}, best_scan={s.best_scan}")
    add("")
    add("Scans")
    for j, (name, n) in enumerate(zip(res.scan_names, res.n_returns), 1):
        st = res.per_scan["status"][:, j - 1]
        used = int(np.sum(res.best_scan == j))
        add(f"  {j}  {name}: {n['leaf']} leaf, {n['wood']} wood, {n['noise']} noise returns; "
            f"chosen for {used} voxels")
        parts = [f"{ladm.SCAN_STATUS[c]} {int(np.sum(st == c))}" for c in sorted(ladm.SCAN_STATUS) if np.any(st == c)]
        add("     voxels by status: " + ", ".join(parts))
        if j <= len(getattr(res, "scan_axes", [])):
            add("     head axis: " + res.scan_axes[j - 1]["text"])
    add("")
    add("Voxels by quality flag")
    for f, name in FLAG_NAMES.items():
        m = res.flag == f
        if np.any(m):
            la = np.nansum(res.lad[m]) * h ** 3
            add(f"  {f} {name:13s} {int(m.sum()):8d} voxels   leaf area {la:10.3f} m2")
    ok = np.isin(res.flag, (FLAG_OK, FLAG_LOW_SAMPLING)) & np.isfinite(res.lad)
    add("")
    add(f"Total leaf area       {res.total_leaf_area():.3f} m2  (flags ok + low_sampling)")
    add(f"                      {res.total_leaf_area((FLAG_OK,)):.3f} m2  (flag ok only)")
    if np.any(ok):
        add(f"Mean LAD              {np.mean(res.lad[ok]):.3f} m2/m3 over {int(ok.sum())} voxels with a value; "
            f"maximum {np.max(res.lad[ok]):.2f}")
    nun = int(np.sum(res.flag == FLAG_UNRESOLVED))
    if nun:
        rows = np.arange(res.lad.size)
        b = np.maximum(res.best_scan - 1, 0)
        nl = res.per_scan["ni_leaf"][rows, b][res.flag == FLAG_UNRESOLVED].sum()
        add(f"Unresolved voxels     {nun}; they hold {int(nl)} leaf returns of the scan with the best view "
            "and are not in the total")
    add("")
    add("Leaf area by height")
    zlow, area, nvox, nnan = vertical_profile(res)
    add("  z_low    z_high   leaf_area_m2  voxels  without_value")
    for z, a, n, nn in zip(zlow, area, nvox, nnan):
        if n:
            add(f"  {z:7.2f}  {z + h:7.2f}  {a:11.3f}  {n:6d}  {nn:6d}")
    if res.warnings:
        add("")
        add("Warnings")
        for w in res.warnings:
            add("  - " + w)
    add("")
    t = res.timings
    add(f"Run time              {t.get('total', 0):.2f} s (counting {t.get('counting', 0):.2f} s)")
    return "\n".join(L) + "\n"


def write_results(res, prefix, per_scan=True, npz=False, matlab_text=False, viewer=True):
    """Write result files; returns the list of paths written."""
    d = os.path.dirname(prefix)
    if d:
        os.makedirs(d, exist_ok=True)
    h = res.voxel_size
    written = []
    centre = res.lo + h / 2.0
    rows = np.arange(res.lad.size)
    b = np.maximum(res.best_scan - 1, 0)
    has = res.best_scan > 0
    pick = lambda k: np.where(has, res.per_scan[k][rows, b], 0)
    cols = {
        "i": res.ijk[:, 0], "j": res.ijk[:, 1], "k": res.ijk[:, 2],
        "x": centre[:, 0], "y": centre[:, 1], "z": centre[:, 2],
        "lad": res.lad, "leaf_area": res.leaf_area, "flag": res.flag.astype(np.int64),
        "best_scan": res.best_scan.astype(np.int64), "percent_explored": res.explored,
        "nt": pick("nt"), "nb": pick("nb"), "ni_leaf": pick("ni_leaf"),
        "ni_wood": pick("ni_wood"), "ni_noise": pick("ni_noise"),
    }
    p = prefix + "_voxels.csv"
    io.write_csv(p, cols, formats={"x": "%.4f", "y": "%.4f", "z": "%.4f"})
    written.append(p)
    if per_scan:
        n, ns = res.per_scan["lad"].shape
        rep = lambda a: np.repeat(a, ns)
        ps = res.per_scan
        cols2 = {"i": rep(res.ijk[:, 0]), "j": rep(res.ijk[:, 1]), "k": rep(res.ijk[:, 2]),
                 "scan": np.tile(np.arange(1, ns + 1), n), "status": ps["status"].ravel().astype(np.int64),
                 "zenith_deg": ps["zenith_deg"].ravel(), "polar_deg": ps["polar_deg"].ravel(),
                 "range": ps["range"].ravel(),
                 "aproj": ps["aproj"].ravel(), "d_bar": ps["d_bar"].ravel(), "G": ps["G"].ravel(),
                 "nt": ps["nt"].ravel(), "nb": ps["nb"].ravel(), "ni_leaf": ps["ni_leaf"].ravel(),
                 "ni_wood": ps["ni_wood"].ravel(), "ni_noise": ps["ni_noise"].ravel(), "ng": ps["ng"].ravel(),
                 "percent_explored": ps["explored"].ravel(), "transmission": ps["transmission"].ravel(),
                 "lad": ps["lad"].ravel()}
        p = prefix + "_scans.csv"
        io.write_csv(p, cols2)
        written.append(p)
    p = prefix + "_summary.txt"
    with open(p, "w") as fh:
        fh.write(summary_text(res))
    written.append(p)
    if npz:
        p = prefix + "_grid.npz"
        np.savez_compressed(p, lad=res.grid("lad"), flag=res.grid("flag", fill=-1),
                            best_scan=res.grid("best_scan", fill=-1), percent_explored=res.grid("explored"),
                            origin=res.origin, voxel_size=h)
        written.append(p)
    if matlab_text:
        p = prefix + "_matlab.txt"
        io.write_matlab_text(p, res.lo, res.per_scan["lad"], res.per_scan["explored"],
                             np.where(res.best_scan > 0, res.best_scan, 1))
        written.append(p)
    if viewer:
        from .viewer import write_viewer
        written.append(write_viewer(res, prefix + "_viewer.html", csv=os.path.basename(prefix) + "_voxels.csv"))
    return written
