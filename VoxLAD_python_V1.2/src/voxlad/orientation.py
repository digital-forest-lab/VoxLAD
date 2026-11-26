"""Rotation axis of a panoramic scanner, found from the returns of the scan.

The pulse count of the "spherical" pattern assumes equal steps in azimuth and
polar angle about the axis around which the scanner head turns.  When the
instrument is not level, that axis is not the z axis of the point cloud.  This
module finds it without any reading from the instrument, from one property that
panoramic scanners share: the returns of one scan line (one sweep of the mirror)
lie on a surface that turns with the head.

Model.  With p the vector from the scanner to a return and rho its length, the
returns of line k satisfy

    (p - D) . m_k = rho * sin(delta) + h

where m_k is a unit vector perpendicular to the head axis a, at head angle
psi_k.  delta = 0 and h = 0 give planes that contain the axis (the meridians of
an ideal scanner); delta is the half-angle by which the lines are cones (1.04
degrees for a RIEGL VZ-400), h a lateral offset of the beam origin, and D an
error of the declared scanner position, perpendicular to the axis.  The axis,
delta, h and D are common to all lines; psi_k is free for each line.  They are
fitted by least squares on the distances of the returns from their line surface.

The lines are taken from the order of the returns in the file: a new line starts
where the direction of consecutive returns jumps.  The method therefore needs a
file that keeps the returns in the order of acquisition, or line by line.
"""
from __future__ import annotations

import numpy as np

_Z = np.array([0.0, 0.0, 1.0])


# ------------------------------------------------------------------ small tools
def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if not n > 0:
        raise ValueError("a direction of zero length was given")
    return v / n


def frame(axis):
    """Right-handed basis (e1, e2, a) with a = axis and e2 as close to +y as possible,
    so that for a vertical axis e1 = +x, e2 = +y."""
    a = unit(axis)
    ref = np.array([0.0, 1.0, 0.0])
    if abs(a @ ref) > 0.999:
        ref = np.array([0.0, 0.0, -1.0]) if a[1] > 0 else np.array([0.0, 0.0, 1.0])
    e2 = unit(ref - (ref @ a) * a)
    return np.cross(e2, a), e2, a


def tilt_and_bearing(axis):
    """Angle of the axis from the vertical (degrees) and the compass-like bearing
    (clockwise from +y, degrees) towards which its upper end leans."""
    a = unit(axis)
    if a[2] < 0:
        a = -a
    tilt = float(np.degrees(np.arccos(min(1.0, a[2]))))
    bearing = float(np.degrees(np.arctan2(a[0], a[1])) % 360.0) if tilt > 1e-9 else 0.0
    return tilt, bearing


def axis_from_tilt(tilt_deg, bearing_deg):
    t, b = np.radians(tilt_deg), np.radians(bearing_deg)
    return np.array([np.sin(t) * np.sin(b), np.sin(t) * np.cos(b), np.cos(t)])


def axis_from_matrix(m):
    """Head axis and scanner position from a registration matrix (scanner frame -> common
    frame; 4 x 4, 3 x 4, or 3 x 3 for the rotation alone).  The scanner's own z axis is its head axis."""
    m = np.asarray(m, dtype=float)
    if m.shape not in ((4, 4), (3, 4), (3, 3)):
        raise ValueError("a registration matrix has 4 x 4, 3 x 4 or 3 x 3 numbers")
    R = m[:3, :3]
    if not np.allclose(R @ R.T, np.eye(3), atol=1e-3) or np.linalg.det(R) < 0:
        raise ValueError("the upper left 3 x 3 block of the registration matrix is not a rotation")
    return unit(R[:, 2]), (m[:3, 3].copy() if m.shape[1] == 4 else None)


def _moved_basis(b1, b2, a0, a):
    """b1, b2 carried along by the smallest rotation that takes a0 to a."""
    k = 1.0 + a0 @ a
    s = a0 + a
    return b1 - (b1 @ a) / k * s, b2 - (b2 @ a) / k * s


# ------------------------------------------------------------------ line model
def _segments(u, jump):
    ch2 = ((u[1:] - u[:-1]) ** 2).sum(axis=1)
    return np.concatenate([[0], np.flatnonzero(ch2 > jump * jump) + 1])


def _moments(F, starts):
    K, f = len(starts), F.shape[1]
    S = np.empty((K, f, f))
    for i in range(f):
        for j in range(i, f):
            S[:, i, j] = S[:, j, i] = np.add.reduceat(F[:, i] * F[:, j], starts)
    return S


class _Fit:
    """Least squares over the lines whose moment matrices are S (K, 5, 5); features (px, py, pz, rho, 1)."""

    names = ("t1", "t2", "delta", "h", "d1", "d2")
    steps = np.array([1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6])

    def __init__(self, S, a0, free_offset=True):
        self.S = S
        self.b1, self.b2, self.a0 = frame(a0)
        self.free = np.array([True, True, True, True, free_offset, free_offset])

    def axis(self, x):
        return unit(self.a0 + x[0] * self.b1 + x[1] * self.b2)

    def W(self, x):
        a = self.axis(x)
        e1, e2 = _moved_basis(self.b1, self.b2, self.a0, a)
        W = np.zeros((3, 5))
        W[0, :3] = e1
        W[0, 4] = -x[4]
        W[1, :3] = e2
        W[1, 4] = -x[5]
        W[2, 3] = -np.sin(x[2])
        W[2, 4] = -x[3]
        return W

    def psi(self, M):
        p, q, r = M[:, 0, 0], M[:, 0, 1], M[:, 1, 1]
        m02, m12, m22 = M[:, 0, 2], M[:, 1, 2], M[:, 2, 2]

        def f(t):
            c, s = np.cos(t), np.sin(t)
            return p * c * c + 2 * q * c * s + r * s * s + 2 * m02 * c + 2 * m12 * s + m22
        t0 = 0.5 * np.arctan2(2 * q, p - r) + np.pi / 2
        t = np.where(f(t0) <= f(t0 + np.pi), t0, t0 + np.pi)
        for _ in range(8):
            c, s = np.cos(t), np.sin(t)
            g = -p * c * s + q * (c * c - s * s) + r * s * c - m02 * s + m12 * c
            H = (r - p) * (c * c - s * s) - 4 * q * c * s - m02 * c - m12 * s
            ok = H > 1e-12 * (p + r)
            t = t - np.where(ok, g / np.where(ok, H, 1.0), 0.0)
        return t, np.maximum(f(t), 0.0)

    def cost(self, x):
        W = self.W(x)
        return self.psi(W @ self.S @ W.T)[1]

    def solve(self, x, iters=40):
        x = np.array(x, dtype=float)
        idx = np.flatnonzero(self.free)
        f0 = float(self.cost(x).sum())
        lam = 1e-3
        H = None
        for _ in range(iters):
            W = self.W(x)
            t, _f = self.psi(W @ self.S @ W.T)
            c = np.stack([np.cos(t), np.sin(t), np.ones_like(t)], axis=1)   # (K, 3)
            dc = np.stack([-np.sin(t), np.cos(t), np.zeros_like(t)], axis=1)
            w = c @ W                                                       # (K, 5) line surface in feature space
            wp = dc @ W
            Sw = (self.S @ w[:, :, None])[:, :, 0]
            Swp = (self.S @ wp[:, :, None])[:, :, 0]
            d = (wp * Swp).sum(axis=1)
            gpsi = (wp * Sw).sum(axis=1)
            dw = []
            for i in idx:
                e = np.zeros(6)
                e[i] = self.steps[i]
                dw.append(c @ ((self.W(x + e) - self.W(x - e)) / (2 * self.steps[i])))   # (K, 5)
            dw = np.stack(dw, axis=1)                                       # (K, m, 5)
            g = (dw @ Sw[:, :, None])[:, :, 0]
            b = (dw @ Swp[:, :, None])[:, :, 0]
            A = dw @ self.S @ dw.transpose(0, 2, 1)
            dd = np.where(d > 0, d, 1.0)
            H = (A - b[:, :, None] * b[:, None, :] / dd[:, None, None]).sum(axis=0)
            G = (g - b * (gpsi / dd)[:, None]).sum(axis=0)
            done = False
            for _try in range(30):
                Hl = H + lam * np.diag(np.diag(H) + 1e-300)
                try:
                    dx = -np.linalg.solve(Hl, G)
                except np.linalg.LinAlgError:
                    lam *= 10
                    continue
                if np.abs(dx).max() < 1e-9:                                       # nothing left to gain
                    done = True
                    break
                x1 = x.copy()
                x1[idx] += dx
                if abs(x1[0]) > 0.3 or abs(x1[1]) > 0.3 or abs(x1[2]) > 0.5:      # stay where the model is sound
                    lam *= 10
                    continue
                f1 = float(self.cost(x1).sum())
                if f1 <= f0:
                    small = (f0 - f1) <= 1e-9 * f0
                    x, f0, lam = x1, f1, max(lam / 5, 1e-9)
                    done = small
                    break
                lam *= 10
            else:
                break
            if done and np.abs(dx).max() < 1e-9:
                break
            # keep the tangent plane centred on the current axis, so that large tilts are followed
            a = self.axis(x)
            e1, e2 = _moved_basis(self.b1, self.b2, self.a0, a)
            self.b1, self.b2, self.a0 = unit(e1), unit(e2), a
            x[0] = x[1] = 0.0
            if done:
                break
        return x, f0, H


# ------------------------------------------------------------------ the estimate
def decimal_quantum(xyz):
    """Rounding step of coordinates written with a fixed number of decimals (0.001 for millimetres);
    0 when no such step is found."""
    n = len(xyz)
    v = np.asarray(xyz[:: max(1, n // 5000)], dtype=float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0
    for k in range(0, 7):
        w = v * 10.0 ** k
        if np.mean(np.abs(w - np.round(w)) < 1e-4) > 0.999:
            return 10.0 ** (-k)
    return 0.0


class AxisEstimate(dict):
    """Result of estimate_axis: a dict with attribute access."""
    __getattr__ = dict.__getitem__


def _sample_blocks(n, max_points, blocks=8):
    if n <= max_points:
        return [(0, n)]
    size = max_points // blocks
    starts = np.linspace(0, n - size, blocks).astype(np.int64)
    return [(int(s), int(s) + size) for s in starts]


def estimate_axis(xyz, position, step=None, max_points=1_200_000, min_range=2.0, jump_deg=2.0,
                  min_line_points=20, start=None):
    """Head axis of a panoramic scanner from its returns, taken in file order.

    xyz        (n, 3) returns in the order of the file
    position   declared scanner position
    step       angular step of the scan (rad), used to judge how well the returns line up
    Returns an AxisEstimate with ``ok`` (bool), ``reason`` (why not, when ok is False),
    ``axis`` (unit vector, upper end), ``tilt_deg``, ``bearing_deg``, ``se_deg`` (standard
    error of the axis), ``cone_deg`` (delta), ``offset`` (vector from the declared position
    to the point where the lines meet, perpendicular to the axis), ``rms_mm`` (distance of
    the returns from their line surface), ``rms_steps`` (the same as a fraction of the
    distance between two lines) and ``lines``, ``points``.
    """
    out = AxisEstimate(ok=False, reason="", axis=_Z.copy(), tilt_deg=0.0, bearing_deg=0.0, se_deg=np.nan,
                       cone_deg=np.nan, offset=np.zeros(3), rms_mm=np.nan, rms_steps=np.nan, lines=0, points=0,
                       median_range=np.nan, step_between=np.nan, step_along=np.nan)
    xyz = np.asarray(xyz)
    pos = np.asarray(position, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 1000:
        out["reason"] = "too few returns"
        return out
    try:
        return _estimate(out, xyz, pos, step, max_points, min_range, jump_deg, min_line_points, start)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        out["ok"] = False
        out["reason"] = "the fit of the scan lines failed"
        return out


def _estimate(out, xyz, pos, step, max_points, min_range, jump_deg, min_line_points, start):
    blocks = _sample_blocks(len(xyz), max_points)
    shift = np.zeros(3)
    best = None
    for _pass in range(4):
        feats, starts, base = [], [], 0
        for a, b in blocks:
            p = np.asarray(xyz[a:b], dtype=float) - (pos + shift)
            rho = np.sqrt((p * p).sum(axis=1))
            keep = rho >= min_range
            p, rho = p[keep], rho[keep]
            if len(p) < 2:
                continue
            st = _segments(p / rho[:, None], np.radians(jump_deg))
            feats.append(np.column_stack([p, rho, np.ones(len(p))]))
            starts.append(st + base)
            base += len(p)
        if not feats:
            out["reason"] = "no return beyond %.1f m from the scanner" % min_range
            return out
        F = np.vstack(feats)
        starts = np.concatenate(starts)
        n = np.diff(np.concatenate([starts, [len(F)]]))
        S = _moments(F, starts)
        mean_rho = S[:, 3, 4] / np.maximum(n, 1)
        sel = n >= min_line_points
        S, n, mean_rho = S[sel], n[sel], mean_rho[sel]
        if len(S) < 30:
            out["reason"] = "fewer than 30 scan lines could be told apart in the file"
            out["lines"] = int(len(S))
            return out
        if best is not None:
            # a later pass, seen from the corrected position: start from the solution in hand
            fit = _Fit(S, best[7])
            x, f, H = fit.solve([0, 0, best[1][2], best[1][3], 0, 0])
        else:
            # starting axes: the vertical, the axis given, and the axis that the free fit of long lines suggests
            cands = [_Z] if start is None else [unit(start), _Z]
            long_ = n >= 200
            if long_.sum() >= 30:
                C = S[long_][:, :3, :3]                     # second moments of p, weighted by rho^2: fine for a start
                nk = np.linalg.eigh(C)[1][:, :, 0]
                cands.append(np.linalg.eigh(np.einsum("ki,kj->ij", nk, nk))[1][:, 0])
            good_enough = (0.25 * (step if step else np.radians(0.04)) * float(np.median(mean_rho))) ** 2 * max(n.sum(), 1)
            fits = []
            for a0 in cands:
                for d0 in (np.radians(0.5), np.radians(3.0)):
                    fit = _Fit(S, a0)
                    x, f, H = fit.solve([0, 0, d0, 0, 0, 0])
                    fits.append((f, x, fit, H))
                    if f < good_enough:                     # the returns are within a quarter of a step of their lines
                        break
                if min(f_[0] for f_ in fits) < good_enough:
                    break
            f, x, fit, H = min(fits, key=lambda q: q[0])
        # leave out lines that fit badly (cut or merged in the wrong place), twice
        for _ in range(2):
            c = fit.cost(x) / n
            good = c <= 16 * np.median(c)
            if good.all() or good.sum() < 30:
                break
            S, n, mean_rho = S[good], n[good], mean_rho[good]
            fit = _Fit(S, fit.axis(x))
            x, f, H = fit.solve([0, 0, x[2], x[3], x[4], x[5]])
        a = fit.axis(x)
        e1, e2 = _moved_basis(fit.b1, fit.b2, fit.a0, a)
        D = x[4] * e1 + x[5] * e2
        best = (f, x, fit, H, S, n, mean_rho, a, D)
        if np.linalg.norm(D) < 0.10:
            break
        shift = shift + D                                   # the lines meet far from the declared position: look again from there
    f, x, fit, H, S, n, mean_rho, a, D = best
    if a[2] < 0:
        a = -a
    dof = max(int(n.sum()) - len(S) - 6, 1)
    rms = float(np.sqrt(f / dof))
    # scatter between interleaved subsets of lines: a standard error that does not rely on the noise being white
    groups = 8
    sub = []
    for gidx in range(groups):
        pick = np.arange(len(S)) % groups == gidx
        if pick.sum() >= 12:
            fg = _Fit(S[pick], fit.axis(x))
            xg, _, _ = fg.solve([0, 0, x[2], x[3], x[4], x[5]], iters=12)
            ag = fg.axis(xg)
            sub.append(ag if ag @ a > 0 else -ag)
    se = np.nan
    if len(sub) >= 4:
        ang = np.arccos(np.clip(np.array(sub) @ a, -1, 1))
        se = float(np.degrees(np.sqrt((ang ** 2).mean() / len(sub))))
    tilt, bearing = tilt_and_bearing(a)
    med_rho = float(np.median(mean_rho))
    out.update(axis=a, tilt_deg=tilt, bearing_deg=bearing, se_deg=se, cone_deg=float(abs(np.degrees(x[2]))),
               offset=shift + D, rms_mm=1000 * rms, lines=int(len(S)), points=int(n.sum()), median_range=med_rho)
    spacing = (step if step else np.radians(0.04)) * med_rho         # distance between two lines at the median range
    q = decimal_quantum(xyz)                                         # rounding of the file adds q / sqrt(12) of white noise
    out["rms_steps"] = float(np.sqrt(max(rms * rms - q * q / 12.0, 0.0)) / spacing)
    if not np.all(np.isfinite(a)) or not np.isfinite(rms):
        out["reason"] = "the fit of the scan lines failed"
    elif not np.isfinite(se):
        out["reason"] = "too few scan lines for an estimate of its uncertainty"
    elif out["rms_steps"] > 0.5:
        out["reason"] = ("the returns do not line up as scan lines seen from the declared scanner position "
                         "(%.1f mm from their line, %.1f of the distance between lines)" % (1000 * rms, out["rms_steps"]))
    elif se > 0.05:
        out["reason"] = "the axis is not determined well enough by this scan (standard error %.2f degrees)" % se
    else:
        out["ok"] = True
        try:
            out["step_between"], out["step_along"] = _measure_steps(
                xyz, pos + out["offset"], a, float(x[2]), float(x[3]), blocks, min_range, step, q)
        except (ValueError, FloatingPointError, IndexError):
            pass
    return out


def _measure_steps(xyz, origin, axis, delta, h, blocks, min_range, step, q, lag=12):
    """Angular spacing of the returns as they stand in the file (rad): between scan lines and along
    a line (between consecutive returns of fully hit surfaces).  nan where it cannot be measured.

    Between lines: each line gets its head angle (median over its returns of the line coordinate of
    the fitted model); the spacing is the median change of that angle over ``lag`` consecutive
    lines, divided by lag.  Twelve lines average over the unequal facets of a polygon mirror, and the
    median is not moved by the few places where a line holds no return.
    """
    if not step:
        return np.nan, np.nan
    e1, e2, a = frame(axis)
    rmin = max(min_range, 2.9 * q / step)                 # the rounding of the file stays below a tenth of a step
    gaps, gaps2 = [], []
    sum_d, sum_n = 0.0, 0
    n_all = n_kept = 0
    for b0, b1 in blocks:
        p = np.asarray(xyz[b0:b1], dtype=float) - origin
        rho = np.sqrt((p * p).sum(axis=1))
        n_all += int((rho >= min_range).sum())
        p, rho = p[rho >= rmin], rho[rho >= rmin]
        n_kept += len(p)
        if len(p) < 200:
            continue
        ca = np.clip((p @ a) / rho, -1.0, 1.0)
        pol = np.arccos(ca)
        sp = np.maximum(np.sqrt(1.0 - ca * ca), 1e-6)
        head = np.unwrap(np.arctan2(p @ e1, p @ e2))
        cone = np.arcsin(np.clip((np.sin(delta) + h / rho) / sp, -1.0, 1.0))
        psi = head + cone                                     # the line coordinate: the same for all returns of a line
        if abs(delta) > 1e-4 or abs(h) > 1e-4:
            alt = head - cone
            if np.median(np.abs(np.diff(alt))) < np.median(np.abs(np.diff(psi))):
                psi = alt
        dpol = np.diff(pol)
        small = np.abs(dpol) < 2.0 * step
        if small.sum() < 100:
            continue
        sg = np.sign(np.median(dpol[small]))
        if sg == 0:
            continue
        # a new line starts where the returns jump back along the line
        span = float(np.percentile(pol, 95) - np.percentile(pol, 5))
        reset = dpol * sg < -max(20.0 * step, 0.25 * span)
        bounds = np.concatenate([[0], np.flatnonzero(reset) + 1, [len(pol)]])
        nline = len(bounds) - 1
        if nline >= lag + 6:
            med = np.array([np.median(psi[bounds[i]:bounds[i + 1]]) for i in range(nline)])
            gaps.append((med[lag:] - med[:-lag]) / lag)
            if nline >= 2 * lag + 3:
                gaps2.append((med[2 * lag:] - med[:-2 * lag]) / (2 * lag))
        # along the lines: runs of ten or more consecutive returns one step apart
        same = np.abs(np.diff(psi)) < 0.4 * step
        adj = same & (dpol * sg > 0.5 * step) & (dpol * sg < 1.5 * step)
        edges = np.flatnonzero(np.diff(np.concatenate([[0], adj.astype(np.int8), [0]])))
        s0, s1 = edges[::2], edges[1::2]
        L = s1 - s0
        big = L >= 10
        if big.any():
            sum_d += float(((pol[s1[big]] - pol[s0[big]]) * sg).sum())
            sum_n += int(L[big].sum())
    if n_kept < 0.6 * n_all:                              # coordinates rounded too coarsely for this
        return np.nan, np.nan
    between = np.nan
    if gaps and gaps2:
        g, g2 = np.concatenate(gaps), np.concatenate(gaps2)
        m, m2 = float(np.median(g)), float(np.median(g2))
        # kept only when the lines advance steadily: lines without a return would stretch the longer lag more
        if g.size >= 60 and g2.size >= 30 and m != 0 and np.median(np.abs(g - m)) < 0.1 * abs(m) \
                and abs(m2 / m - 1.0) < 0.01:
            between = abs(m)
    along = sum_d / sum_n if sum_n >= 2000 else np.nan
    return between, along


# ------------------------------------------------------------------ the axis of a scan
def resolve_axis(scan, matlab=False):
    """Decide which head axis the model uses for a scan, and check the scan against its description.

    Returns a dict: ``axis`` (unit vector, or None when the axis is the vertical), ``tilt_deg``,
    ``bearing_deg``, ``source`` ("vertical", "given", "matrix" or "estimated"), ``text`` (one
    line for the summary), ``estimate`` (the AxisEstimate from the scan lines, when the returns
    were analysed) and ``warnings``.  The result is kept on the scan, so the returns are analysed once.
    """
    spec = getattr(scan, "axis", "auto")
    key = (spec if isinstance(spec, str) else tuple(np.round(np.asarray(spec, dtype=float), 12)), bool(matlab))
    cached = getattr(scan, "_axis_info", None)
    if cached is not None and cached.get("_key") == key:
        return cached
    info = {"axis": None, "tilt_deg": 0.0, "bearing_deg": 0.0, "source": "vertical", "text": "", "estimate": None,
            "warnings": [], "_key": key}
    name = getattr(scan, "name", "") or "scan"
    matrix = getattr(scan, "matrix", None)

    def use(a, source):
        a = unit(a)
        if a[2] < 0:
            a = -a
        info["tilt_deg"], info["bearing_deg"] = tilt_and_bearing(a)
        info["axis"] = None if (a[0] == 0.0 and a[1] == 0.0) else a
        info["source"] = source

    def lean(a=None):
        if a is None and info["axis"] is None:
            return "vertical"
        t, b = (info["tilt_deg"], info["bearing_deg"]) if a is None else tilt_and_bearing(a)
        if t < 0.005:
            return "vertical to within 0.01 deg"
        return "tilted %.2f deg from the vertical, towards bearing %.0f deg" % (t, b)

    if matlab:
        info["text"] = "taken as vertical (MATLAB mode)"
        scan._axis_info = info
        return info
    # the scan lines tell the axis, and also whether the declared position and angular step fit the returns
    est = estimate_axis(scan.xyz, scan.position, scan.step_az)
    info["estimate"] = est
    found = " (found from %s scan lines, +-%.3f deg)" % (format(est.lines, ","), max(est.se_deg, 0.001)) if est.ok else ""
    if not isinstance(spec, str):
        use(spec, "given")
        info["text"] = lean() + " (as given)"
    elif spec == "vertical":
        info["text"] = "taken as vertical (as set)"
    elif spec == "matrix":
        use(axis_from_matrix(matrix)[0], "matrix")
        info["text"] = lean() + " (from the registration matrix)"
    elif scan.pattern != "spherical":
        if matrix is not None:
            use(axis_from_matrix(matrix)[0], "matrix")
            info["text"] = lean() + " (from the registration matrix; used for the limits of the scan only)"
        else:
            info["text"] = "taken as vertical; the pulse count of the uniform pattern does not depend on it"
            if est.ok:
                info["text"] += " (the scan lines turn about an axis %.1f deg from the vertical)" % est.tilt_deg
    elif est.ok:
        use(est.axis, "estimated")
        info["text"] = lean() + found
    elif matrix is not None:
        use(axis_from_matrix(matrix)[0], "matrix")
        info["text"] = lean() + " (from the registration matrix; not found from the returns: %s)" % est.reason
    else:
        info["text"] = "taken as vertical; the tilt could not be found from the returns: " + est.reason
    if est.ok:
        if info["source"] in ("given", "matrix") and scan.pattern == "spherical":
            used = _Z if info["axis"] is None else info["axis"]
            diff = float(np.degrees(np.arccos(min(1.0, abs(float(used @ est.axis))))))
            if diff > 0.2:
                info["warnings"].append(
                    "scan '%s': the head axis %s is %.2f deg from the axis found from the scan lines (%s)."
                    % (name, "given" if info["source"] == "given" else "of the registration matrix", diff, lean(est.axis)))
        elif info["source"] == "estimated" and matrix is not None:
            diff = float(np.degrees(np.arccos(min(1.0, abs(float(axis_from_matrix(matrix)[0] @ est.axis))))))
            if diff > 0.2:
                info["warnings"].append(
                    "scan '%s': the head axis found from the scan lines is %.2f deg from the z axis of its "
                    "registration matrix; the axis found from the scan lines was used." % (name, diff))
        off = float(np.linalg.norm(est.offset))
        if off > 0.10:
            info["warnings"].append(
                "scan '%s': its scan lines meet %.2f m from the declared scanner position (dx, dy, dz = "
                "%.2f, %.2f, %.2f m, measured across the axis the lines turn about). Check 'position'."
                % (name, off, est.offset[0], est.offset[1], est.offset[2]))
        rb = est.step_between / scan.step_az if np.isfinite(est.step_between) else np.nan
        ra = est.step_along / scan.step_zen if np.isfinite(est.step_along) else np.nan
        info["step_ratio"] = (float(rb), float(ra))
        if (np.isfinite(rb) and abs(rb - 1) > 0.02) or (np.isfinite(ra) and abs(ra - 1) > 0.02):
            parts = []
            if np.isfinite(rb):
                parts.append("%.5g deg between scan lines" % np.degrees(est.step_between))
            if np.isfinite(ra):
                parts.append("%.5g deg along them" % np.degrees(est.step_along))
            dens = 1.0 / ((rb if np.isfinite(rb) else 1.0) * (ra if np.isfinite(ra) else 1.0))
            info["warnings"].append(
                "scan '%s': the returns are spaced %s, where the declared angular step is %.5g x %.5g deg. "
                "With the spacing found, N_T would be %.0f%% %s. Check the angular step."
                % (name, " and ".join(parts), np.degrees(scan.step_az), np.degrees(scan.step_zen),
                   abs(100 * (dens - 1)), "lower" if dens < 1 else "higher"))
    scan._axis_info = info
    return info


def tilt_effect(tilt_deg, elevation_deg):
    """Largest relative error of N_T (fraction) at a given elevation above the scanner's horizon
    when a head axis tilted by ``tilt_deg`` is taken as vertical."""
    th, t = np.radians(90.0 - elevation_deg), np.radians(tilt_deg)
    return float(max(abs(np.sin(th - t) / np.sin(th) - 1), abs(np.sin(th + t) / np.sin(th) - 1)))
