"""Configuration file (TOML) -> Settings and scans."""
from __future__ import annotations

import os

import numpy as np

from . import io
from .io import InputError
from .model import Scan, Settings

try:
    import tomllib as _toml            # Python 3.11+
except ImportError:  # pragma: no cover
    try:
        import tomli as _toml
    except ImportError:
        _toml = None

_GRID_KEYS = {"voxel_size", "bounds", "z_range", "voxels"}
_MODEL_KEYS = {"leaf_angle", "beam_correction", "jensen_fit", "unresolved", "best_scan",
               "min_percent_explored", "min_range_voxels", "compat", "engine", "threads"}
_CLASS_KEYS = {"leaf", "wood", "noise", "unlisted"}
_SCAN_KEYS = {"file", "name", "format", "position", "pattern", "angular_step_deg", "angular_step_rad",
              "left_edge", "right_edge", "zenith_range_deg", "columns", "skip_rows", "delimiter",
              "returns", "class_field", "axis", "matrix", "matrix_file"}
_OUTPUT_KEYS = {"prefix", "per_scan", "npz", "matlab_text", "viewer", "open_viewer"}
_TOP_KEYS = {"grid", "model", "classes", "scanner", "scan", "output", "advanced"}

TEMPLATE = '''\
# VoxLAD configuration file (TOML).  Lines starting with # are comments.
# Lengths are in metres.  All scans must be in the same coordinate frame, z up.

[grid]
voxel_size = 0.3                # side of the cubic voxels
# bounds = [-5, 16, -1, 4, 28, 7]   # xmin ymin zmin xmax ymax zmax; default: extent of the leaf returns
# z_range = [-1.0, 12.0]        # only evaluate voxels between these heights
# voxels = "leaf"               # "leaf": voxels holding a leaf return; "all": every voxel of the grid

[model]
leaf_angle = "spherical"        # erectophile, planophile, plagiophile, spherical, uniform,
                                # or a list of fractions for equal classes from 0 to 90 degrees
# beam_correction = 1.0         # multiplies the leaf returns in a voxel (beam-size correction H)
# min_percent_explored = 15     # voxels sampled below this are flagged low_sampling
# jensen_fit = "paper"          # "paper" (Eq. 7), "exact", "matlab", "matlab-scaled", "none"
# unresolved = "nan"            # "nan", "beer" or "zero" for voxels with no real solution
# best_scan = "entering"        # "entering": most pulses reaching the voxel; "fraction": highest share
# compat = "matlab"             # reproduce the MATLAB version (see the manual)

[classes]                       # class codes used in your files
leaf = [1]
wood = [2]
noise = [3]
# unlisted = "error"            # other codes: "error", "noise" or "drop"

[scanner]                       # applies to every scan; any key may be repeated inside a [[scan]]
pattern = "spherical"           # "spherical": constant azimuth and zenith steps (panoramic scanners)
                                # "uniform": same spacing in both directions (window scanners, ILRIS-3D)
angular_step_deg = 0.04         # one number, or [azimuth, zenith]; or use angular_step_rad
zenith_range_deg = [30, 130]    # vertical field of view, measured from straight up
# axis = "auto"                 # the axis the scanner head turns about.  "auto": found from the scan lines,
                                # so that an instrument that was not level is handled; "vertical"; or [x, y, z]

[[scan]]
file = "scan1.txt"
position = [0.0, 0.0, 0.0]      # scanner optical centre in the common frame
# columns = ["x", "y", "z", "class"]   # for text files; add "return_number" if present
# left_edge = [-1.0, 0.0]       # x y of a point on the left limit of the scan, seen from the scanner
# right_edge = [1.0, 0.0]       # x y of a point on the right limit; omit both for 360 degrees
# matrix_file = "scan1.dat"     # registration matrix of the scan (4 x 4 numbers, scanner frame to common frame):
                                # gives the axis when it cannot be found from the returns, and the position if absent

[[scan]]
file = "scan2.txt"
position = [10.0, 0.0, 0.0]

[output]
prefix = "voxlad"               # files are written as <prefix>_voxels.csv, <prefix>_summary.txt, ...
# per_scan = true               # also write <prefix>_scans.csv (one line per voxel and scan)
# npz = false                   # also write <prefix>_grid.npz (3-D arrays for Python)
# matlab_text = false           # also write <prefix>_matlab.txt in the layout of the MATLAB version
# viewer = true                 # write <prefix>_viewer.html, an interactive 3-D view of the result
# open_viewer = true            # open it in the web browser at the end of a command-line run
'''


def _check_keys(table, allowed, where):
    extra = set(table) - allowed
    if extra:
        raise InputError(f"unknown setting(s) {sorted(extra)} in {where}. Allowed: {sorted(allowed)}")


def _step(table, where):
    has_deg = "angular_step_deg" in table
    has_rad = "angular_step_rad" in table
    if has_deg and has_rad:
        raise InputError(f"{where}: give angular_step_deg or angular_step_rad, not both")
    if not (has_deg or has_rad):
        return None
    v = np.atleast_1d(np.asarray(table["angular_step_deg" if has_deg else "angular_step_rad"], dtype=float))
    if v.size == 1:
        v = np.repeat(v, 2)
    if v.size != 2 or np.any(v <= 0):
        raise InputError(f"{where}: the angular step must be one positive number or [azimuth, zenith]")
    return v * (np.pi / 180.0) if has_deg else v


def load_config(path):
    """Read a configuration file.  Returns (settings, scan_specs, output dict).

    scan_specs are plain dicts; call load_scans() to read the point clouds.
    """
    if _toml is None:
        raise InputError("reading TOML needs Python 3.11 or the 'tomli' package (pip install tomli)")
    try:
        with open(path, "rb") as fh:
            cfg = _toml.load(fh)
    except FileNotFoundError:
        raise InputError(f"configuration file not found: {path}") from None
    except _toml.TOMLDecodeError as err:
        raise InputError(f"{path}: not valid TOML: {err}") from None
    _check_keys(cfg, _TOP_KEYS, "the configuration file")
    grid = cfg.get("grid", {})
    model = cfg.get("model", {})
    classes = cfg.get("classes", {})
    scanner = cfg.get("scanner", {})
    output = cfg.get("output", {})
    _check_keys(grid, _GRID_KEYS, "[grid]")
    _check_keys(model, _MODEL_KEYS, "[model]")
    _check_keys(classes, _CLASS_KEYS, "[classes]")
    _check_keys(scanner, _SCAN_KEYS - {"file", "name", "position"}, "[scanner]")
    _check_keys(output, _OUTPUT_KEYS, "[output]")
    if "voxel_size" not in grid:
        raise InputError("[grid] voxel_size is required")
    settings = Settings(advanced=dict(cfg.get("advanced", {})), **{**grid, **model})
    settings.resolved()            # validates
    classes = {"leaf": classes.get("leaf", [1]), "wood": classes.get("wood", [2]),
               "noise": classes.get("noise", [3]), "unlisted": classes.get("unlisted", "error")}
    scans = cfg.get("scan", [])
    if not scans:
        raise InputError("no [[scan]] entry in the configuration file")
    base = os.path.dirname(os.path.abspath(path))
    specs = []
    for n, sc in enumerate(scans, 1):
        where = f"[[scan]] number {n}"
        _check_keys(sc, _SCAN_KEYS, where)
        if "file" not in sc:
            raise InputError(f"{where}: 'file' is required")
        spec = {**scanner, **sc}
        spec["path"] = sc["file"] if os.path.isabs(sc["file"]) else os.path.join(base, sc["file"])
        spec["name"] = sc.get("name") or os.path.splitext(os.path.basename(sc["file"]))[0]
        if "matrix_file" in spec and not os.path.isabs(spec["matrix_file"]):
            spec["matrix_file"] = os.path.join(base, spec["matrix_file"])
        spec["_where"] = where
        spec["_classes"] = classes
        specs.append(spec)
    out = {"prefix": output.get("prefix", "voxlad"), "per_scan": output.get("per_scan", True),
           "npz": output.get("npz", False), "matlab_text": output.get("matlab_text", False),
           "viewer": output.get("viewer", True), "open_viewer": output.get("open_viewer", True)}
    if not os.path.isabs(out["prefix"]):
        out["prefix"] = os.path.join(base, out["prefix"])
    return settings, specs, out


def _matrix(spec, where):
    """Registration matrix of a scan from 'matrix' (numbers) or 'matrix_file' (a text file of numbers)."""
    if "matrix" in spec and "matrix_file" in spec:
        raise InputError(f"{where}: give matrix or matrix_file, not both")
    if "matrix" in spec:
        try:
            v = np.asarray(spec["matrix"], dtype=float).ravel()
        except (ValueError, TypeError):
            raise InputError(f"{where}: matrix must hold numbers only") from None
        src = "matrix"
    elif "matrix_file" in spec:
        f = spec["matrix_file"]
        if not os.path.isfile(f):
            raise InputError(f"{where}: matrix file not found: {f}")
        tokens = open(f, errors="replace").read().replace(",", " ").replace(";", " ").split()
        nums = []
        for t in tokens:
            try:
                nums.append(float(t))
            except ValueError:
                pass
        v = np.asarray(nums)
        src = os.path.basename(f)
    else:
        return None
    if v.size not in (9, 12, 16):
        raise InputError(f"{where}: {src} holds {v.size} numbers; a registration matrix has 16 (4 x 4), "
                         "12 (3 x 4) or 9 (3 x 3, rotation only)")
    return v.reshape((3, 3) if v.size == 9 else (-1, 4))


def load_scan(spec):
    """Read one scan described by a spec dict from load_config()."""
    path, where = spec["path"], spec.get("_where", "scan")
    if not os.path.isfile(path):
        raise InputError(f"{where}: file not found: {path}")
    fmt = spec.get("format", "auto")
    if fmt == "auto":
        fmt = io.detect_format(path)
    header = {}
    columns = spec.get("columns", ["x", "y", "z", "class"])
    if fmt == "voxlad-asc":
        header, data = io.read_voxlad_asc(path)
        cols = {"x": data[:, 0], "y": data[:, 1], "z": data[:, 2], "class": data[:, 3]}
    elif fmt == "text":
        cols = io.read_text(path, columns, spec.get("skip_rows", 0), spec.get("delimiter"))
    elif fmt == "las":
        cols = io.read_las(path, spec.get("class_field", "classification"))
    elif fmt == "npy":
        cols = io.read_npy(path, columns)
    else:
        raise InputError(f"{where}: unknown format '{fmt}' (auto, voxlad-asc, text, las, npy)")
    xyz = np.column_stack([cols["x"], cols["y"], cols["z"]]).astype(float)
    n_file = xyz.shape[0]
    keep = np.ones(n_file, dtype=bool)
    returns = spec.get("returns", "auto")
    if returns not in ("auto", "first", "all"):
        raise InputError(f"{where}: returns must be 'auto', 'first' or 'all'")
    if "return_number" in cols and returns in ("auto", "first"):
        keep &= np.asarray(cols["return_number"]) <= 1
    elif returns == "first":
        raise InputError(f"{where}: returns = \"first\" needs a return_number column")
    cls, keep_cls = io.map_classes(cols["class"], spec["_classes"], spec["_classes"]["unlisted"],
                                   where=f"{os.path.basename(path)}: ")
    keep &= keep_cls
    # scanner description: the configuration overrides the file header
    matrix = _matrix(spec, where)
    position = spec.get("position", header.get("position"))
    if position is None and matrix is not None and matrix.shape[1] == 4:
        position = matrix[:3, 3]
    if position is None:
        raise InputError(f"{where}: 'position' (scanner x y z) is required")
    axis = spec.get("axis", "auto")
    step = _step(spec, where)
    if step is None:
        if "step_az" not in header:
            raise InputError(f"{where}: angular_step_deg or angular_step_rad is required")
        step = (header["step_az"], header["step_zen"])
    pattern = spec.get("pattern")
    if pattern is None:
        raise InputError(f"{where}: 'pattern' is required (\"spherical\" or \"uniform\"); "
                         "it can be set once under [scanner]")
    if "left_edge" in spec or "right_edge" in spec:
        left, right = spec.get("left_edge"), spec.get("right_edge")
    else:
        left, right = header.get("left_edge"), header.get("right_edge")
    zr = spec.get("zenith_range_deg", [0.0, 180.0])
    if len(zr) != 2 or not (0 <= zr[0] < zr[1] <= 180):
        raise InputError(f"{where}: zenith_range_deg must be [min, max] within 0..180")
    try:
        scan = Scan(xyz[keep], cls[keep], position, pattern, float(step[0]), float(step[1]),
                    left, right, np.radians(zr[0]), np.radians(zr[1]), spec["name"], axis=axis, matrix=matrix)
    except InputError as err:
        raise InputError(f"{where}: {err}") from None
    scan.n_in_file = n_file
    scan.path = path
    scan.format = fmt
    scan.npts_header = header.get("npts_header")
    return scan


def load_scans(specs, progress=None):
    scans = []
    for spec in specs:
        sc = load_scan(spec)
        if progress:
            progress(f"read {sc.xyz.shape[0]} returns from {os.path.basename(sc.path)} ({sc.format})")
        scans.append(sc)
    return scans
