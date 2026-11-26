"""Reading point clouds and writing results."""
from __future__ import annotations

import os

import numpy as np

LEAF, WOOD, NOISE = 1, 2, 3


class InputError(ValueError):
    """A problem with an input file or setting, reported without a traceback."""


def _loadtxt(fh, **kw):
    try:
        return np.loadtxt(fh, ndmin=2, **kw)
    except ValueError as err:
        raise InputError(f"could not read numeric columns: {err}") from None


def read_voxlad_asc(path):
    """Read the original VoxLAD text format.

    line 1  x y z of the scanner
    line 2  x y of a point on the left edge of the scan (seen from the scanner)
    line 3  x y of a point on the right edge; give the scanner x y on both
            lines for a 360-degree scan
    line 4  'd' or 'r' (degrees or radians) and the angular step
    line 5  number of returns (informative; the file is read to its end)
    then    x y z class, one return per line

    Returns (header dict, data array (n, 4)).
    """
    with open(path) as fh:
        try:
            pos = [float(v) for v in fh.readline().split()]
            left = [float(v) for v in fh.readline().split()]
            right = [float(v) for v in fh.readline().split()]
            unit, step = fh.readline().split()
            npts = int(float(fh.readline().split()[0]))
            if len(pos) != 3 or len(left) != 2 or len(right) != 2 or unit not in ("d", "r"):
                raise ValueError
        except (ValueError, IndexError):
            raise InputError(
                f"{path}: not a VoxLAD .asc file (expected a 5-line header: scanner "
                "x y z / left edge x y / right edge x y / 'd' or 'r' and angular step / "
                "number of returns)") from None
        data = _loadtxt(fh)
    if data.shape[1] < 4:
        raise InputError(f"{path}: expected at least 4 columns (x y z class), found {data.shape[1]}")
    step = float(step) * (np.pi / 180.0 if unit == "d" else 1.0)
    header = {"position": pos, "left_edge": left, "right_edge": right,
              "step_az": step, "step_zen": step, "npts_header": npts}
    return header, data


def looks_like_voxlad_asc(path):
    try:
        with open(path) as fh:
            lines = [fh.readline().split() for _ in range(5)]
        return ([len(t) for t in lines[:3]] == [3, 2, 2] and len(lines[3]) == 2
                and lines[3][0] in ("d", "r") and len(lines[4]) == 1)
    except (OSError, UnicodeDecodeError):
        return False


def read_text(path, columns, skip_rows=0, delimiter=None):
    """Plain text or CSV: one return per line; ``columns`` names the columns in
    order ('x', 'y', 'z', 'class', optionally 'return_number'; use '-' or any
    other name for columns to ignore)."""
    if delimiter is None:
        with open(path) as fh:
            for _ in range(skip_rows):
                fh.readline()
            first = fh.readline()
        delimiter = "," if "," in first else (";" if ";" in first else None)
    data = _loadtxt(path, skiprows=skip_rows, delimiter=delimiter, comments="#")
    return _pick_columns(data, columns, path)


def _pick_columns(data, columns, path):
    names = [c.lower() for c in columns]
    for need in ("x", "y", "z", "class"):
        if need not in names:
            raise InputError(f"{path}: the column list {columns} has no '{need}'")
    if data.shape[1] < len(names):
        raise InputError(f"{path}: {len(names)} columns named but the file has {data.shape[1]}")
    out = {k: data[:, names.index(k)] for k in ("x", "y", "z", "class")}
    if "return_number" in names:
        out["return_number"] = data[:, names.index("return_number")]
    return out


def read_las(path, class_field="classification"):
    """LAS or LAZ through the optional laspy package."""
    try:
        import laspy
    except ImportError:
        raise InputError("reading LAS/LAZ files needs the 'laspy' package: pip install laspy "
                         "(and 'lazrs' for .laz)") from None
    las = laspy.read(path)
    try:
        cls = np.asarray(las[class_field])
    except Exception:
        raise InputError(f"{path}: no point field named '{class_field}'") from None
    out = {"x": np.asarray(las.x, dtype=float), "y": np.asarray(las.y, dtype=float),
           "z": np.asarray(las.z, dtype=float), "class": cls}
    try:
        out["return_number"] = np.asarray(las.return_number)
    except Exception:
        pass
    try:
        # a file sorted in space loses the order of acquisition, which the search for the tilt of the
        # scanner relies on: the time stamps bring it back (the order of the returns changes nothing else)
        t = np.asarray(las.gps_time, dtype=float)
        if t.size == out["x"].size and np.ptp(t) > 0 and np.any(np.diff(t) < 0):
            order = np.argsort(t, kind="stable")
            out = {k: v[order] for k, v in out.items()}
    except Exception:
        pass
    return out


def read_npy(path, columns):
    """NumPy .npy holding an (n, k) array, or .npz with arrays x, y, z, class."""
    if str(path).lower().endswith(".npz"):
        z = np.load(path)
        missing = [k for k in ("x", "y", "z", "class") if k not in z]
        if missing:
            raise InputError(f"{path}: arrays {missing} not found")
        out = {k: np.asarray(z[k]) for k in ("x", "y", "z", "class")}
        if "return_number" in z:
            out["return_number"] = np.asarray(z["return_number"])
        return out
    arr = np.load(path)
    if arr.ndim != 2:
        raise InputError(f"{path}: expected a 2-D array")
    return _pick_columns(arr, columns, path)


def detect_format(path):
    ext = os.path.splitext(str(path))[1].lower()
    if ext in (".las", ".laz"):
        return "las"
    if ext in (".npy", ".npz"):
        return "npy"
    return "voxlad-asc" if looks_like_voxlad_asc(path) else "text"


def map_classes(raw, classes, unlisted="error", where=""):
    """Translate the class codes found in a file to LEAF / WOOD / NOISE.

    Returns (codes int8, keep mask).  Codes that are not listed are an error
    by default, or become noise ('noise'), or are removed ('drop').
    """
    raw = np.asarray(raw)
    rawi = np.rint(raw.astype(float)).astype(np.int64)
    out = np.zeros(raw.shape, dtype=np.int8)
    for code, key in ((LEAF, "leaf"), (WOOD, "wood"), (NOISE, "noise")):
        vals = np.atleast_1d(np.asarray(classes.get(key, []), dtype=np.int64))
        if vals.size:
            out[np.isin(rawi, vals)] = code
    keep = np.ones(raw.shape, dtype=bool)
    unknown = out == 0
    if np.any(unknown):
        found = np.unique(rawi[unknown])
        if unlisted == "error":
            raise InputError(
                f"{where}class codes {found.tolist()[:20]} are not assigned to leaf, wood or noise. "
                "List them under [classes], or set unlisted = \"noise\" or \"drop\".")
        if unlisted == "noise":
            out[unknown] = NOISE
        elif unlisted == "drop":
            keep = ~unknown
        else:
            raise InputError("classes.unlisted must be 'error', 'noise' or 'drop'")
    return out, keep


def _fmt_matlab(v, fmt):
    if np.isnan(v):
        return "NaN"
    if np.isinf(v):
        return "-Inf" if v < 0 else "Inf"
    return fmt % v


def write_matlab_text(path, lo, lad, explored, best):
    """Result file laid out as the MATLAB script writes it: per voxel
    'x_min y_min z_min' then 'LAD percent_explored' for each scan, then the
    number of the best scan."""
    with open(path, "w") as fh:
        fh.write("x_min y_min z_min LAD_v_scan_1 percent_explored_scan1 ...nscan best_scan\n")
        for i in range(lo.shape[0]):
            parts = ["%.3f %.3f %.3f " % tuple(lo[i])]
            for s in range(lad.shape[1]):
                parts.append("%s %s " % (_fmt_matlab(lad[i, s], "%.2f"), _fmt_matlab(explored[i, s], "%.0f")))
            parts.append("%d\n" % best[i])
            fh.write("".join(parts))


def write_csv(path, columns, header_comment=None, formats=None):
    """columns: dict name -> 1-D array.  Integer columns are written as
    integers, others with 6 significant digits unless ``formats`` says otherwise."""
    names = list(columns)
    formats = formats or {}
    fmts = []
    for k in names:
        a = np.asarray(columns[k])
        fmts.append(formats.get(k, "%d" if np.issubdtype(a.dtype, np.integer) else "%.6g"))
    arr = np.column_stack([np.asarray(columns[k], dtype=float) for k in names])
    with open(path, "w") as fh:
        if header_comment:
            for line in header_comment.splitlines():
                fh.write("# " + line + "\n")
        fh.write(",".join(names) + "\n")
        np.savetxt(fh, arr, fmt=fmts, delimiter=",")
