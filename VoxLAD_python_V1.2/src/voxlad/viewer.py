"""Interactive 3-D viewer of the result.

The viewer is a single HTML file that opens in any current web browser: the
voxel values are stored inside the file and drawn with WebGL, so nothing has to
be installed and the file can be sent to a colleague as it is.  Voxels can be
rotated, moved and zoomed, coloured by leaf area density, percent explored,
quality flag or scan used, filtered by value and cut along x, y or z; pointing
at a voxel shows its values.
"""
from __future__ import annotations

import base64
import html
import json
import os
import sys
import tempfile

import numpy as np

from . import __version__
from ._viewer_template import TEMPLATE
from .io import InputError


def _pack(a, dtype):
    a = np.ascontiguousarray(a, dtype=np.dtype(dtype).newbyteorder("<"))
    return {"dtype": np.dtype(dtype).name, "b64": base64.b64encode(a.tobytes()).decode("ascii")}


def _payload(title, voxel_size, origin, shape, ijk, lad, flag, scan, explored, nt, nb, nil, niw, nin,
             scans=(), csv=None, file=None):
    ijk = np.asarray(ijk)
    if ijk.size and (ijk.min() < 0 or ijk.max() > 65535):
        raise InputError("the viewer handles grids of up to 65,535 voxels along an axis")
    clip = lambda a: np.clip(np.nan_to_num(np.asarray(a, dtype=float)), -2 ** 31 + 1, 2 ** 31 - 1)
    return {
        "title": str(title), "version": __version__, "voxel_size": float(voxel_size),
        "origin": [float(v) for v in origin], "shape": [int(v) for v in shape], "n": int(len(lad)),
        "scans": [{"name": str(n), "position": [float(v) for v in p]} for n, p in scans],
        "csv": csv, "file": file,
        "arrays": {
            "ijk": _pack(ijk.reshape(-1, 3), np.uint16), "lad": _pack(lad, np.float32),
            "explored": _pack(explored, np.float32), "flag": _pack(np.clip(flag, 0, 4), np.uint8),
            "scan": _pack(np.clip(scan, 0, 255), np.uint8), "nt": _pack(clip(nt), np.int32),
            "nb": _pack(clip(nb), np.int32), "nil": _pack(clip(nil), np.int32),
            "niw": _pack(clip(niw), np.int32), "nin": _pack(clip(nin), np.int32),
        },
    }


def _payload_from_result(res, title=None, csv=None, file=None):
    rows = np.arange(res.lad.size)
    b = np.maximum(res.best_scan - 1, 0)
    has = res.best_scan > 0
    pick = lambda k: np.where(has, res.per_scan[k][rows, b], 0)
    positions = getattr(res, "scan_positions", None) or []
    scans = list(zip(res.scan_names, positions)) if len(positions) == len(res.scan_names) else []
    return _payload(title or "VoxLAD result", res.voxel_size, res.origin, res.shape, res.ijk, res.lad, res.flag,
                    res.best_scan, res.explored, pick("nt"), pick("nb"), pick("ni_leaf"), pick("ni_wood"),
                    pick("ni_noise"), scans, csv, file)


def _payload_from_csv(path, title=None):
    """Rebuild what the viewer needs from a <prefix>_voxels.csv file."""
    try:
        t = np.genfromtxt(path, delimiter=",", names=True)
    except (OSError, ValueError) as err:
        raise InputError(f"cannot read {path}: {err}") from None
    t = np.atleast_1d(t)
    need = ("i", "j", "k", "x", "y", "z", "lad", "flag", "best_scan", "percent_explored")
    names = t.dtype.names or ()
    missing = [c for c in need if c not in names]
    if missing:
        raise InputError(f"{path} is not a VoxLAD voxel table: column(s) {missing} not found")
    ijk = np.column_stack([t["i"], t["j"], t["k"]]).astype(np.int64)
    xyz = np.column_stack([t["x"], t["y"], t["z"]])
    h = None
    for a in range(3):                                   # voxel size from two voxels that differ along one axis
        d = ijk[:, a] - ijk[0, a]
        m = d != 0
        if np.any(m):
            h = float(np.median((xyz[m, a] - xyz[0, a]) / d[m]))
            break
    if h is None or not h > 0:
        raise InputError(f"{path}: the voxel size cannot be worked out from a single voxel")
    origin = xyz[0] - (ijk[0] + 0.5) * h
    col = lambda c: t[c] if c in names else np.zeros(t.size)
    name = os.path.basename(path)
    stem = name[:-len("_voxels.csv")] if name.endswith("_voxels.csv") else os.path.splitext(name)[0]
    return _payload(title or stem, round(h, 6), np.round(origin, 6), ijk.max(axis=0) + 1, ijk,
                    t["lad"], t["flag"].astype(int), t["best_scan"].astype(int), t["percent_explored"],
                    col("nt"), col("nb"), col("ni_leaf"), col("ni_wood"), col("ni_noise"), csv=name, file=stem)


def build_html(payload):
    data = json.dumps(payload, separators=(",", ":")).replace("<", "\\u003c")
    return TEMPLATE.replace("__TITLE__", html.escape(payload["title"])).replace("__DATA__", data)


def write_viewer(res, path, title=None, csv=None):
    """Write the viewer of a Result to ``path`` (an .html file).  Returns the path."""
    stem = os.path.basename(path)
    stem = stem[:-len("_viewer.html")] if stem.endswith("_viewer.html") else os.path.splitext(stem)[0]
    payload = _payload_from_result(res, title or stem, csv, stem)
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(build_html(payload))
    return path


def write_viewer_from_csv(csv_path, path=None, title=None):
    """Write a viewer from a <prefix>_voxels.csv file written by an earlier run."""
    if path is None:
        base = csv_path[:-len("_voxels.csv")] if csv_path.endswith("_voxels.csv") else os.path.splitext(csv_path)[0]
        path = base + "_viewer.html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(build_html(_payload_from_csv(csv_path, title)))
    return path


def can_open_browser():
    """False on a machine without a graphical session (a remote shell, a cluster node)."""
    if sys.platform.startswith(("win", "darwin", "cygwin")):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def open_in_browser(path):
    """Open an HTML file in the default web browser.  Returns True if a browser was started."""
    if not can_open_browser():
        return False
    import webbrowser
    from pathlib import Path
    try:
        return bool(webbrowser.open(Path(path).resolve().as_uri()))
    except Exception:
        return False


def view(res, path=None, title=None, open=True):
    """Write the viewer of a Result and open it in the web browser.

    Without ``path`` the file goes to the temporary folder.  Returns the path.
    """
    if path is None:
        fd, path = tempfile.mkstemp(prefix="voxlad_", suffix="_viewer.html")
        os.close(fd)
    write_viewer(res, path, title)
    if open:
        open_in_browser(path)
    return path
