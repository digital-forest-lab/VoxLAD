"""The VoxLAD pipeline: scans in, leaf area density per voxel out."""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from . import geometry as geo
from . import lad as ladm
from . import leafangle as la
from . import legacy, orientation, traverse
from .io import LEAF, WOOD, NOISE, InputError

# quality flag of the final value of a voxel
FLAG_OK = 0
FLAG_LOW_SAMPLING = 1    # fewer than min_percent_explored of the pulses reached the voxel
FLAG_UNRESOLVED = 2      # no scan gives a real solution (voxel too dense for the sampling)
FLAG_OCCLUDED = 3        # covered by a scan but no pulse reached it
FLAG_NOT_SEEN = 4        # outside the field of view (or too close) for every scan
FLAG_NAMES = {FLAG_OK: "ok", FLAG_LOW_SAMPLING: "low_sampling", FLAG_UNRESOLVED: "unresolved",
              FLAG_OCCLUDED: "occluded", FLAG_NOT_SEEN: "not_seen"}

MAX_VOXELS = 400_000_000


@dataclass
class Scan:
    """One scan: returns, scanner position and sampling pattern.

    xyz        (n, 3) coordinates of the returns, metres, z up, all scans in
               the same frame
    cls        (n,) class of each return: 1 leaf, 2 wood, 3 noise
    position   scanner optical centre (x, y, z) in the same frame
    pattern    "spherical" (constant azimuth and zenith steps: panoramic
               scanners) or "uniform" (same spacing in both directions
               everywhere: window scanners such as the ILRIS-3D)
    step_az, step_zen   angular steps in radians
    left_edge, right_edge   optional (x, y) points on the left and right limits
               of the scan as seen from the scanner; None for 360 degrees
    zenith_min, zenith_max  limits of the vertical field of view, radians from
               the zenith (0 = straight up); for a tilted scanner, from its axis
    axis       the axis around which the scanner head turns: "auto" (found from
               the scan lines when the returns are in the order of the file,
               else vertical), "vertical", "matrix" (z axis of ``matrix``), or a
               vector (x, y, z) in the common frame
    matrix     optional registration matrix of the scan (scanner frame to common
               frame, 4 x 4); gives the axis when it cannot be found from the returns
    """
    xyz: np.ndarray
    cls: np.ndarray
    position: Sequence[float]
    pattern: str
    step_az: float
    step_zen: Optional[float] = None
    left_edge: Optional[Sequence[float]] = None
    right_edge: Optional[Sequence[float]] = None
    zenith_min: float = 0.0
    zenith_max: float = float(np.pi)
    name: str = ""
    axis: object = "auto"
    matrix: Optional[Sequence] = None

    def __post_init__(self):
        self.xyz = np.ascontiguousarray(self.xyz, dtype=np.float64).reshape(-1, 3)
        self.cls = np.ascontiguousarray(self.cls, dtype=np.int8).ravel()
        if self.cls.shape[0] != self.xyz.shape[0]:
            raise InputError("xyz and cls must have the same number of returns")
        self.position = np.asarray(self.position, dtype=float).reshape(3)
        if self.step_zen is None:
            self.step_zen = self.step_az
        if self.pattern not in ("spherical", "uniform"):
            raise InputError("scan pattern must be 'spherical' or 'uniform'")
        if not (self.step_az > 0 and self.step_zen > 0):
            raise InputError("angular steps must be positive")
        if not np.all(np.isfinite(self.xyz)):
            raise InputError(f"scan {self.name!r}: non-finite coordinates")
        if (self.left_edge is None) != (self.right_edge is None):
            raise InputError("give both left_edge and right_edge, or neither")
        if self.matrix is not None:
            try:
                self.matrix = np.asarray(self.matrix, dtype=float)
                orientation.axis_from_matrix(self.matrix)
            except (ValueError, TypeError) as err:
                raise InputError(f"scan {self.name!r}: {err}") from None
        if isinstance(self.axis, str):
            if self.axis not in ("auto", "vertical", "matrix"):
                raise InputError('axis must be "auto", "vertical", "matrix" or a vector [x, y, z]')
            if self.axis == "matrix" and self.matrix is None:
                raise InputError(f'scan {self.name!r}: axis = "matrix" needs a registration matrix')
        else:
            try:
                a = np.asarray(self.axis, dtype=float).ravel()
                if a.size != 3 or not np.all(np.isfinite(a)):
                    raise ValueError
                self.axis = orientation.unit(a)
            except (ValueError, TypeError):
                raise InputError('axis must be "auto", "vertical", "matrix" or a vector [x, y, z]') from None


@dataclass
class Settings:
    voxel_size: float = 0.3
    bounds: Optional[Sequence[float]] = None      # xmin ymin zmin xmax ymax zmax
    z_range: Optional[Sequence[float]] = None     # zmin zmax
    voxels: str = "leaf"                          # "leaf" or "all"
    leaf_angle: object = "spherical"              # name or list of class fractions
    beam_correction: object = None                # None, a number, or one factor per class
    jensen_fit: Optional[str] = None              # paper | matlab | matlab-scaled | none
    unresolved: Optional[str] = None              # nan | zero | beer
    best_scan: str = "entering"                   # entering | fraction
    min_percent_explored: float = 15.0
    min_range_voxels: float = 2.0
    compat: str = "none"                          # none | matlab
    engine: Optional[str] = None                  # traversal | matlab
    threads: Optional[int] = None
    advanced: dict = field(default_factory=dict)

    def resolved(self):
        """Fill in the settings that depend on ``compat``."""
        if self.compat not in ("none", "matlab"):
            raise InputError("compat must be 'none' or 'matlab'")
        m = self.compat == "matlab"
        opt = {
            "engine": self.engine or ("matlab" if m else "traversal"),
            "jensen_fit": self.jensen_fit or ("matlab" if m else "paper"),
            "unresolved": self.unresolved or ("zero" if m else "nan"),
            "face_rule": "matlab" if m else "exact",
            "coverage": "matlab" if m else "exact",
            "g_folded": not m,
            "lad_rule": "matlab" if m else "new",
            "grid_rule": "matlab" if m else "new",
            "filter_points_by_height": m,
        }
        unknown = set(self.advanced) - set(opt)
        if unknown:
            raise InputError(f"unknown advanced option(s): {sorted(unknown)}")
        opt.update(self.advanced)
        if opt["engine"] not in ("traversal", "matlab"):
            raise InputError("engine must be 'traversal' or 'matlab'")
        if opt["jensen_fit"] not in ladm.JENSEN_FITS:
            raise InputError(f"jensen_fit must be one of {ladm.JENSEN_FITS}")
        if opt["unresolved"] not in ("nan", "zero", "beer"):
            raise InputError("unresolved must be 'nan', 'zero' or 'beer'")
        if opt["jensen_fit"] == "exact" and opt["lad_rule"] == "matlab":
            raise InputError('jensen_fit = "exact" cannot be combined with compat = "matlab"')
        if self.best_scan not in ("entering", "fraction"):
            raise InputError("best_scan must be 'entering' or 'fraction'")
        if self.voxels not in ("leaf", "all"):
            raise InputError("voxels must be 'leaf' or 'all'")
        if not self.voxel_size > 0:
            raise InputError("voxel_size must be positive")
        return opt


@dataclass
class Result:
    settings: Settings
    options: dict
    scan_names: list
    voxel_size: float
    origin: np.ndarray            # lower corner of voxel (0, 0, 0)
    shape: tuple                  # voxels along x, y, z
    ijk: np.ndarray               # (n, 3) grid indices of the evaluated voxels
    lo: np.ndarray                # (n, 3) their lower corners
    per_scan: dict                # name -> (n, nscan) array
    lad: np.ndarray               # (n,) final leaf area density, m2/m3
    flag: np.ndarray              # (n,) quality flag
    best_scan: np.ndarray         # (n,) 1-based scan number used, 0 if none
    explored: np.ndarray          # (n,) percent explored for that scan
    g_table: np.ndarray
    beam_factor: float
    n_returns: list               # per scan: dict of counts by class
    scan_positions: list = field(default_factory=list)    # scanner x y z, per scan
    scan_axes: list = field(default_factory=list)         # per scan: head axis used and how it was obtained
    timings: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    @property
    def leaf_area(self):
        """Leaf area per voxel (m2): LAD times voxel volume."""
        return self.lad * self.voxel_size ** 3

    def total_leaf_area(self, flags=(FLAG_OK, FLAG_LOW_SAMPLING)):
        m = np.isin(self.flag, flags) & np.isfinite(self.lad)
        return float(np.sum(self.lad[m]) * self.voxel_size ** 3)

    def grid(self, what="lad", fill=np.nan):
        """Dense (nx, ny, nz) array of a final quantity ('lad', 'flag',
        'best_scan' or 'explored'); voxels that were not evaluated get ``fill``."""
        src = getattr(self, what)
        out = np.full(self.shape, fill, dtype=float if what in ("lad", "explored") else np.int16)
        out[self.ijk[:, 0], self.ijk[:, 1], self.ijk[:, 2]] = src
        return out


def _beam_factor(spec, g):
    if spec is None:
        return 1.0
    if isinstance(spec, str):
        if spec.lower() in ("ilris3d-20m", "ilris"):
            spec = la.H_ILRIS3D_20M
        else:
            raise InputError("beam_correction must be a number, a list, or 'ilris3d-20m'")
    arr = np.atleast_1d(np.asarray(spec, dtype=float))
    if arr.size == 1:
        return float(arr[0])
    if arr.size != g.size:
        raise InputError(f"beam_correction has {arr.size} values but the leaf angle "
                         f"distribution has {g.size} classes")
    # sum_i g_i*H_i, accumulated in the same order as the MATLAB loop
    tot = 0.0
    for gi, hi in zip(g, arr):
        tot = tot + gi * hi
    return float(tot)


def _matlab_axis(b, bmax, h):
    n = int(np.floor((bmax - b) / h + 1e-9)) + 1
    return b + np.arange(n) * h


def _build_grid(scans, s, opt):
    """Returns origin, shape, and per-axis lower-bound arrays."""
    h = s.voxel_size
    if s.bounds is not None:
        b = np.asarray(s.bounds, dtype=float).ravel()
        if b.size != 6 or np.any(b[3:] <= b[:3]):
            raise InputError("bounds must be [xmin, ymin, zmin, xmax, ymax, zmax] with max > min")
        lo, hi = b[:3].copy(), b[3:].copy()
    else:
        if opt["grid_rule"] == "matlab":
            pts = [sc.xyz for sc in scans if sc.xyz.shape[0]]
        else:
            pts = [sc.xyz[sc.cls == LEAF] for sc in scans if np.any(sc.cls == LEAF)]
            if s.z_range is not None:
                pts = [p[(p[:, 2] >= s.z_range[0]) & (p[:, 2] <= s.z_range[1])] for p in pts]
                pts = [p for p in pts if p.shape[0]]
        if not pts:
            raise InputError("no leaf returns found: nothing to voxelise (check [classes] and z_range)")
        allp = np.vstack(pts)
        lo, hi = np.floor(allp.min(axis=0)), np.ceil(allp.max(axis=0))
        if s.z_range is not None and opt["grid_rule"] != "matlab":
            lo[2] = max(lo[2], np.floor(s.z_range[0])) if np.isfinite(s.z_range[0]) else lo[2]
            hi[2] = min(hi[2], np.ceil(s.z_range[1])) if np.isfinite(s.z_range[1]) else hi[2]
    if opt["grid_rule"] == "matlab" and s.bounds is None:
        axes = [_matlab_axis(lo[a], hi[a], h) for a in range(3)]
    else:
        axes = [lo[a] + np.arange(max(1, int(np.ceil((hi[a] - lo[a]) / h - 1e-9)))) * h for a in range(3)]
    shape = tuple(len(a) for a in axes)
    nvox = int(np.prod([float(v) for v in shape]))
    if nvox > MAX_VOXELS:
        raise InputError(
            f"the grid would have {shape[0]} x {shape[1]} x {shape[2]} = {nvox:.3g} voxels. "
            "Set [grid] bounds to the volume of interest or use a larger voxel_size.")
    return lo, shape, axes


def compute(scans, settings=None, progress=None):
    """Run VoxLAD on a list of Scan objects.  Returns a Result."""
    s = settings or Settings()
    opt = s.resolved()
    say = progress or (lambda msg: None)
    if not scans:
        raise InputError("no scans given")
    h = float(s.voxel_size)
    t_start = time.perf_counter()
    timings = {}
    warns = []
    g = la.leaf_angle_table(s.leaf_angle)
    if opt["g_folded"]:
        tot = g.sum()
        if not tot > 0:
            raise InputError("leaf angle class fractions sum to zero")
        if abs(tot - 1.0) > 1e-3:
            warns.append(f"leaf angle class fractions summed to {tot:.4f}; they were rescaled to 1")
        if abs(tot - 1.0) > 1e-6:
            g = g / tot
    beam = _beam_factor(s.beam_correction, g)

    n_returns = [{"leaf": int(np.sum(sc.cls == LEAF)), "wood": int(np.sum(sc.cls == WOOD)),
                  "noise": int(np.sum(sc.cls == NOISE))} for sc in scans]
    work = scans
    if opt["filter_points_by_height"] and s.z_range is not None:
        # MATLAB drops returns outside (min_canopy_height, max_canopy_height) at reading
        work = []
        for sc in scans:
            k = (sc.xyz[:, 2] < s.z_range[1]) & (sc.xyz[:, 2] > s.z_range[0])
            work.append(Scan(sc.xyz[k], sc.cls[k], sc.position, sc.pattern, sc.step_az, sc.step_zen,
                             sc.left_edge, sc.right_edge, sc.zenith_min, sc.zenith_max, sc.name))

    origin, shape, axes = _build_grid(work, s, opt)
    nx, ny, nz = shape
    nvox = nx * ny * nz
    say(f"grid {nx} x {ny} x {nz} voxels of {h:g} m, origin {origin.tolist()}")

    # --- which voxels to evaluate -------------------------------------------
    def index_of(xyz):
        if opt["grid_rule"] == "matlab":
            ijk = np.stack([np.searchsorted(axes[a], xyz[:, a], side="right") - 1 for a in range(3)], axis=1)
            inside = np.all((ijk >= 0) & (ijk < np.array(shape)), axis=1)
            # a return past the last lower bound belongs to the last voxel only if within it
            for a in range(3):
                inside &= xyz[:, a] < axes[a][-1] + h
            return np.where(inside, (ijk[:, 0] * ny + ijk[:, 1]) * nz + ijk[:, 2], -1)
        return traverse.point_voxel_index(xyz, origin, h, shape)

    if s.voxels == "all":
        sel_flat = np.arange(nvox, dtype=np.int64)
    else:
        occ = np.zeros(nvox, dtype=bool)
        for sc in work:
            leaf = sc.xyz[sc.cls == LEAF]
            if s.z_range is not None and s.bounds is None and opt["grid_rule"] != "matlab":
                leaf = leaf[(leaf[:, 2] >= s.z_range[0]) & (leaf[:, 2] <= s.z_range[1])]
            f = index_of(leaf)
            occ[f[f >= 0]] = True
        sel_flat = np.flatnonzero(occ)
    nsel = sel_flat.size
    if nsel == 0:
        raise InputError("no voxel contains a leaf return; check class codes, bounds and z_range")
    lookup = np.full(nvox, -1, dtype=np.int32)
    lookup[sel_flat] = np.arange(nsel, dtype=np.int32)
    ijk = np.stack(np.unravel_index(sel_flat, shape), axis=1).astype(np.int32)
    lo = np.stack([axes[a][ijk[:, a]] for a in range(3)], axis=1)
    centre = lo + h / 2.0
    say(f"{nsel} voxels to evaluate, {len(work)} scan(s), engine '{opt['engine']}'")

    nscan = len(work)
    f8 = lambda: np.full((nsel, nscan), np.nan)
    ps = {k: f8() for k in ("zenith_deg", "polar_deg", "range", "aproj", "d_bar", "G", "transmission", "explored", "lad")}
    for k in ("nt", "nb", "ni_leaf", "ni_wood", "ni_noise", "ng"):
        ps[k] = np.zeros((nsel, nscan), dtype=np.int64)
    ps["status"] = np.zeros((nsel, nscan), dtype=np.int8)
    entering = np.zeros((nsel, nscan))

    timings["setup"] = time.perf_counter() - t_start
    # --- the axis around which each scanner head turned -----------------------
    t0 = time.perf_counter()
    scan_axes = []
    for sc in work:
        info = orientation.resolve_axis(sc, matlab=opt["coverage"] == "matlab")
        scan_axes.append(info)
        warns.extend(info["warnings"])
        say(f"scan '{sc.name}': head axis {info['text']}")
    timings["axis"] = time.perf_counter() - t0
    t_count = t_geom = 0.0
    for j, sc in enumerate(work):
        t0 = time.perf_counter()
        ax = scan_axes[j]["axis"]
        if ax is not None and scan_axes[j]["tilt_deg"] > 20.0 and sc.left_edge is not None \
                and (len(sc.left_edge) < 3 or len(sc.right_edge) < 3):
            raise InputError(
                f"scan '{sc.name}': its head axis is {scan_axes[j]['tilt_deg']:.0f} deg from the vertical, so the limits "
                "of the scan cannot be given by two coordinates. Give left_edge and right_edge as [x, y, z] points.")
        # ---- counts ----------------------------------------------------------
        if opt["engine"] == "matlab":
            cnt = legacy.count_matlab(sc.xyz, sc.cls, lo, h, sc.position)
            nb, ni = cnt[:, 0], cnt[:, 1:4]
        else:
            nb = traverse.count_blocked(sc.xyz, sc.position, origin, h, shape, lookup, nsel, threads=s.threads)
            ni = traverse.count_hits(sc.xyz, sc.cls, origin, h, shape, lookup, nsel)
        t1 = time.perf_counter()
        t_count += t1 - t0
        # ---- geometry --------------------------------------------------------
        aproj, nfaces = geo.projected_area(lo, h, sc.position, face_rule=opt["face_rule"])
        theta = geo.view_zenith(centre, sc.position)
        v = centre - sc.position
        rng = np.sqrt((v ** 2).sum(axis=1))
        with np.errstate(divide="ignore", invalid="ignore"):
            dbar = h ** 3 / aproj
        nt = geo.pulse_count(aproj, centre, sc.position, sc.pattern, sc.step_az, sc.step_zen, axis=ax)
        if opt["coverage"] == "matlab":
            cov = geo.covered_matlab(centre, h, sc.position, sc.left_edge, sc.right_edge,
                                     sc.pattern, sc.zenith_min)
            near = np.zeros(nsel, dtype=bool)
        else:
            cov = geo.covered(lo, h, sc.position, sc.left_edge, sc.right_edge, sc.zenith_min, sc.zenith_max,
                              axis=ax)
            near = (rng < s.min_range_voxels * h) | (nfaces == 0)
            if sc.pattern == "spherical":
                near |= geo.axis_distance(v, ax) < h                  # voxel on the scan's pole
        G = la.g_function(theta, g, folded=opt["g_folded"])
        delta = None if opt["jensen_fit"] == "exact" else ladm.jensen_delta(h, dbar, opt["jensen_fit"])
        nil, niw, nin = ni[:, 0], ni[:, 1], ni[:, 2]
        if opt["lad_rule"] == "matlab":
            # scans MATLAB skips keep their zero-initialised values
            z = lambda a: np.where(cov, a, 0)
            nt_, nb_, nil_, niw_, nin_ = z(nt), z(nb), z(nil), z(niw), z(nin)
            ng = nt_ - nb_ - nil_ - niw_ - nin_
            lad_s = ladm.invert_lad_matlab(z(dbar), ladm.jensen_delta(h, z(dbar), opt["jensen_fit"]),
                                           np.where(cov, G, 0), ng, nil_ * beam)
            with np.errstate(divide="ignore", invalid="ignore"):
                expl = (nt_ - nb_) / nt_ * 100.0
                tr = ng / (ng + nil_ * beam)
            status = np.where(cov, ladm.OK, ladm.NOT_COVERED).astype(np.int8)
            ent = (nt_ - nb_).astype(float)
            nt, nb, nil, niw, nin = nt_, nb_, nil_, niw_, nin_
            Gout = np.where(cov, np.real(G), np.nan)
        else:
            ng = nt - nb - nil - niw - nin
            ent = nt - nb
            status = np.full(nsel, ladm.OK, dtype=np.int8)
            status[ent <= 0] = ladm.OCCLUDED
            status[nt < 1] = ladm.NO_PULSES
            status[near] = ladm.TOO_CLOSE
            status[~cov] = ladm.NOT_COVERED
            usable = status == ladm.OK
            if opt["jensen_fit"] == "exact":
                with np.errstate(divide="ignore", invalid="ignore"):
                    u_abs = np.abs(v) / rng[:, None]
                lad_s, st = ladm.invert_lad_exact(u_abs, h, G, ng, nil * beam, d_bar=dbar)
                if opt["unresolved"] == "zero":
                    lad_s = np.where(st == ladm.NO_GAP, 0.0, lad_s)
            else:
                lad_s, st = ladm.invert_lad(dbar, delta, G, ng, nil * beam, unresolved=opt["unresolved"])
            status = np.where(usable, st, status).astype(np.int8)
            lad_s = np.where(usable, lad_s, np.nan)
            with np.errstate(divide="ignore", invalid="ignore"):
                expl = np.where(nt >= 1, ent / nt * 100.0, np.nan)
                tr = np.where(usable & (ng > 0), ng / (ng + nil * beam), np.where(usable, 0.0, np.nan))
            ent = np.where(usable, ent, 0.0)
            Gout = G
        ps["status"][:, j] = status
        ps["zenith_deg"][:, j] = np.degrees(theta)
        ps["polar_deg"][:, j] = np.degrees(geo.axis_angles(v, ax)[0])
        ps["range"][:, j] = rng
        ps["aproj"][:, j] = aproj
        ps["d_bar"][:, j] = dbar
        ps["G"][:, j] = Gout
        ps["nt"][:, j] = np.where(np.isfinite(nt), nt, 0).astype(np.int64)
        ps["nb"][:, j] = nb
        ps["ni_leaf"][:, j] = nil
        ps["ni_wood"][:, j] = niw
        ps["ni_noise"][:, j] = nin
        ps["ng"][:, j] = np.where(np.isfinite(ng), ng, 0).astype(np.int64)
        ps["explored"][:, j] = expl
        ps["transmission"][:, j] = tr
        ps["lad"][:, j] = lad_s
        entering[:, j] = ent
        t_geom += time.perf_counter() - t1
        say(f"scan {j + 1}/{nscan} '{sc.name}': {sc.xyz.shape[0]} returns, counted in {t1 - t0:.2f} s")
    timings["counting"] = t_count
    timings["geometry_and_inversion"] = t_geom

    # --- choose a scan per voxel --------------------------------------------
    rows = np.arange(nsel)
    status = ps["status"]
    if opt["lad_rule"] == "matlab":
        best = np.argmax(entering, axis=1)                     # first maximum, as MATLAB max()
        lad = ps["lad"][rows, best]
        expl = ps["explored"][rows, best]
        flag = np.where(np.isfinite(lad), FLAG_OK, FLAG_NOT_SEEN).astype(np.int8)
        best1 = best + 1
    else:
        solved = status == ladm.OK
        entered = np.isin(status, (ladm.OK, ladm.NO_GAP, ladm.NO_SOLUTION))
        if s.best_scan == "fraction":
            score = np.where(entered, ps["explored"], -np.inf)
        else:
            score = np.where(entered, entering, -np.inf)
        # prefer a scan that gives a solution; otherwise the best-sampled scan that
        # reached the voxel (its value is nan, or the fallback asked for with `unresolved`)
        score_solved = np.where(solved, score, -np.inf)
        any_solved = solved.any(axis=1)
        best = np.where(any_solved, np.argmax(score_solved, axis=1), np.argmax(score, axis=1))
        any_entered = entered.any(axis=1)
        lad = np.where(any_entered, ps["lad"][rows, best], np.nan)
        expl = np.where(any_entered, ps["explored"][rows, best], np.nan)
        seen = np.isin(status, (ladm.OK, ladm.NO_GAP, ladm.NO_SOLUTION, ladm.OCCLUDED)).any(axis=1)
        flag = np.full(nsel, FLAG_NOT_SEEN, dtype=np.int8)
        flag[seen] = FLAG_OCCLUDED
        flag[any_entered] = FLAG_UNRESOLVED
        ok = any_entered & (status[rows, best] == ladm.OK)
        flag[ok] = FLAG_OK
        flag[ok & (expl < s.min_percent_explored)] = FLAG_LOW_SAMPLING
        best1 = np.where(any_entered, best + 1, 0)

        # sanity check of the declared angular step
        with np.errstate(divide="ignore", invalid="ignore"):
            over = (ps["nb"] + ps["ni_leaf"] + ps["ni_wood"] + ps["ni_noise"]) > 1.2 * ps["nt"] + 3
        cov_n = np.isin(status, (ladm.OK, ladm.NO_GAP, ladm.NO_SOLUTION, ladm.OCCLUDED))
        for j, sc in enumerate(work):
            ncov = int(cov_n[:, j].sum())
            if ncov >= 20:
                frac = float((over[:, j] & cov_n[:, j]).sum()) / ncov
                if frac > 0.10:
                    warns.append(
                        f"scan '{sc.name}': in {100 * frac:.0f}% of voxels more returns were found than "
                        "pulses the declared angular step allows. Check angular_step, the scan pattern, "
                        "and that the file holds one return per pulse.")
    timings["total"] = time.perf_counter() - t_start
    return Result(settings=s, options=opt, scan_names=[sc.name or f"scan{j + 1}" for j, sc in enumerate(work)],
                  voxel_size=h, origin=np.asarray(origin, dtype=float), shape=shape, ijk=ijk, lo=lo,
                  per_scan=ps, lad=lad, flag=flag, best_scan=best1.astype(np.int32), explored=expl,
                  g_table=g, beam_factor=beam, n_returns=n_returns, timings=timings, warnings=warns,
                  scan_positions=[np.asarray(sc.position, dtype=float) for sc in work], scan_axes=scan_axes)
