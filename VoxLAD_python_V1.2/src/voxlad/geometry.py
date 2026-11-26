"""Scanner-voxel geometry: projected area, theoretical pulse count, coverage.

Everything here is vectorised over voxels for one scan.  Voxels are axis
aligned cubes given by their lower corner ``lo`` (n, 3) and side ``h``.
"""
from __future__ import annotations

import numpy as np

# Vertex numbering of the MATLAB code (get_vox_vertex.m), 0-based here:
#   0:(0,0,0) 1:(0,0,1) 2:(1,0,1) 3:(1,0,0) 4:(0,1,0) 5:(0,1,1) 6:(1,1,1) 7:(1,1,0)
VERTEX_OFFSETS = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0],
                           [0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]], dtype=float)
# faces as cyclic vertex lists: axis -> (low-side face, high-side face)
FACES = {0: ((0, 1, 5, 4), (3, 2, 6, 7)),     # x = min, x = max
         1: ((0, 1, 2, 3), (4, 5, 6, 7)),     # y = min, y = max
         2: ((0, 4, 7, 3), (1, 5, 6, 2))}     # z = min, z = max


def voxel_vertices(lo, h):
    """(n, 8, 3) vertex coordinates; the upper faces are lo + h, as in MATLAB."""
    lo = np.asarray(lo, dtype=float)
    hi = lo + h
    v = np.empty((lo.shape[0], 8, 3))
    for k in range(8):
        for a in range(3):
            v[:, k, a] = hi[:, a] if VERTEX_OFFSETS[k, a] else lo[:, a]
    return v


def _quad_area(p):
    """Area of planar quads p (n, 4, 3): four triangles about the centroid."""
    c = p.mean(axis=1)
    area = np.zeros(p.shape[0])
    for i in range(4):
        cp = np.cross(c - p[:, i], c - p[:, (i + 1) % 4])
        area += 0.5 * np.sqrt((cp ** 2).sum(axis=1))
    return area


def projected_area(lo, h, origin, face_rule="exact"):
    """Cross-section A_proj of each voxel seen from ``origin`` (paper, Fig. 3).

    The faces turned towards the scanner are projected along the lines of
    sight onto the plane through the voxel centre perpendicular to the
    scanner-centre line, and their areas are summed.

    face_rule = "exact"   a face counts only if the scanner is in front of it
                          (one, two or three faces)
              = "matlab"  always three faces, the nearer one of each pair, as in
                          the MATLAB code; slightly too large when the scanner
                          lies between the planes of two opposite faces
    Returns (aproj, n_faces).
    """
    lo = np.asarray(lo, dtype=float)
    o = np.asarray(origin, dtype=float)
    n = lo.shape[0]
    hi = lo + h
    centre = lo + h / 2.0
    verts = voxel_vertices(lo, h)
    v = centre - o                              # scanner -> centre
    d2 = (v ** 2).sum(axis=1)
    aproj = np.zeros(n)
    nfaces = np.zeros(n, dtype=np.int8)
    with np.errstate(divide="ignore", invalid="ignore"):
        for a in range(3):
            near_low = np.abs(o[a] - lo[:, a]) < np.abs(o[a] - hi[:, a])
            if face_rule == "matlab":
                use = np.ones(n, dtype=bool)
            else:
                use = (o[a] < lo[:, a]) | (o[a] > hi[:, a])
            idx = np.where(near_low[:, None], np.array(FACES[a][0])[None, :],
                           np.array(FACES[a][1])[None, :])
            fv = np.take_along_axis(verts, idx[:, :, None].repeat(3, axis=2), axis=1)
            w = fv - o                                           # scanner -> vertex
            t = d2[:, None] / (w * v[:, None, :]).sum(axis=2)
            proj = o + w * t[:, :, None]
            area = _quad_area(proj)
            aproj += np.where(use, area, 0.0)
            nfaces += use.astype(np.int8)
    return aproj, nfaces


def matlab_round(x):
    """MATLAB round(): halves away from zero (numpy rounds halves to even)."""
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def is_vertical(axis):
    """True when ``axis`` is None or the z axis: the level-instrument formulas apply as they are."""
    if axis is None:
        return True
    a = np.asarray(axis, dtype=float)
    return bool(a[0] == 0.0 and a[1] == 0.0 and a[2] != 0.0)


def axis_frame(axis):
    """Right-handed basis (e1, e2, a) of a scanner whose head turns about ``axis``.
    For a vertical axis e1 = +x and e2 = +y, so that bearings keep their meaning."""
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    ref = np.array([0.0, 1.0, 0.0])
    if abs(a @ ref) > 0.999:
        ref = np.array([0.0, 0.0, -1.0]) if a[1] > 0 else np.array([0.0, 0.0, 1.0])
    e2 = ref - (ref @ a) * a
    e2 = e2 / np.linalg.norm(e2)
    return np.cross(e2, a), e2, a


def axis_distance(v, axis=None):
    """Distance of the points ``v`` (vectors from the scanner) from the line of the head axis."""
    v = np.asarray(v, dtype=float)
    if is_vertical(axis):
        return np.sqrt(v[..., 0] ** 2 + v[..., 1] ** 2)
    a = axis_frame(axis)[2]
    along = v @ a
    return np.sqrt(np.maximum((v ** 2).sum(axis=-1) - along ** 2, 0.0))


def axis_angles(v, axis=None):
    """Polar angle from the head axis and bearing about it (rad) of the vectors ``v`` from the scanner.
    With a vertical axis these are the zenith angle and the bearing clockwise from +y."""
    v = np.asarray(v, dtype=float)
    d = np.sqrt((v ** 2).sum(axis=-1))
    with np.errstate(divide="ignore", invalid="ignore"):
        if is_vertical(axis):
            return np.arccos(np.clip(v[..., 2] / d, -1.0, 1.0)), np.arctan2(v[..., 0], v[..., 1])
        e1, e2, a = axis_frame(axis)
        return np.arccos(np.clip((v @ a) / d, -1.0, 1.0)), np.arctan2(v @ e1, v @ e2)


def pulse_count(aproj, centre, origin, pattern, step_az, step_zen, axis=None):
    """Theoretical number of pulses N_T aimed at each voxel.

    pattern = "spherical": constant steps in azimuth and in polar angle about
        the axis around which the scanner head turns, so the spacing across the
        beam at range d is 2*d*tan(step_zen/2) along a scan line and
        2*d_axis*tan(step_az/2) between lines, d_axis being the distance of
        the voxel centre from the line of the axis.  ``axis`` is that axis as a
        vector; None means vertical, d_axis then being the horizontal range.
    pattern = "uniform": the same spacing 2*d*tan(step/2) in both directions
        whatever the direction (window scanners such as the ILRIS-3D).
    """
    v = np.asarray(centre, dtype=float) - np.asarray(origin, dtype=float)
    d = np.sqrt((v ** 2).sum(axis=1))
    dxy = axis_distance(v, axis)
    l_lat = 2.0 * d * np.tan(step_zen / 2.0)
    if pattern == "uniform":
        l_lon = 2.0 * d * np.tan(step_az / 2.0)
    elif pattern == "spherical":
        l_lon = 2.0 * dxy * np.tan(step_az / 2.0)
    else:
        raise ValueError("pattern must be 'spherical' or 'uniform'")
    with np.errstate(divide="ignore", invalid="ignore"):
        return matlab_round(aproj / (l_lat * l_lon))


def view_zenith(centre, origin):
    """Zenith angle (rad) of the direction scanner -> voxel centre."""
    v = np.asarray(centre, dtype=float) - np.asarray(origin, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.arccos(v[:, 2] / np.sqrt((v ** 2).sum(axis=1)))


def _bearing(dx, dy):
    """Clockwise angle from +y (rad), i.e. a compass bearing if y is north."""
    return np.arctan2(dx, dy)


def in_sector(x, y, origin, left_edge, right_edge):
    """True where (x, y) lies in the sector swept clockwise (seen from above)
    from the left edge direction to the right edge direction."""
    if left_edge is None or right_edge is None:
        return np.ones(np.shape(x), dtype=bool)
    ox, oy = float(origin[0]), float(origin[1])
    lx, ly = left_edge[0] - ox, left_edge[1] - oy
    rx, ry = right_edge[0] - ox, right_edge[1] - oy
    if (lx == 0 and ly == 0) or (rx == 0 and ry == 0):
        return np.ones(np.shape(x), dtype=bool)        # edge = scanner: 360 degrees
    return _between(_bearing(np.asarray(x) - ox, np.asarray(y) - oy), _bearing(lx, ly), _bearing(rx, ry))


def _between(bearing, b_left, b_right):
    width = np.mod(b_right - b_left, 2 * np.pi)
    if width == 0:
        return np.ones(np.shape(bearing), dtype=bool)
    return np.mod(bearing - b_left, 2 * np.pi) <= width


def in_sector_about(v, origin, left_edge, right_edge, axis=None):
    """As in_sector, for vectors ``v`` from the scanner, with bearings taken about the head axis.
    The edges are the half-planes that hold the axis and the horizontal direction from the
    scanner to each edge point."""
    v = np.asarray(v, dtype=float)
    if left_edge is None or right_edge is None:
        return np.ones(v.shape[:-1], dtype=bool)
    o = np.asarray(origin, dtype=float)
    el = np.array([left_edge[0] - o[0], left_edge[1] - o[1], left_edge[2] - o[2] if len(left_edge) > 2 else 0.0])
    er = np.array([right_edge[0] - o[0], right_edge[1] - o[1], right_edge[2] - o[2] if len(right_edge) > 2 else 0.0])
    if not el.any() or not er.any():
        return np.ones(v.shape[:-1], dtype=bool)
    return _between(axis_angles(v, axis)[1], axis_angles(el, axis)[1], axis_angles(er, axis)[1])


def covered(lo, h, origin, left_edge=None, right_edge=None,
            zenith_min=0.0, zenith_max=np.pi, axis=None):
    """Voxels wholly inside the scan's azimuth sector and zenith range.

    A voxel that straddles a limit is not used, because only part of its
    theoretical pulses were emitted.  All eight vertices are tested.  The
    limits belong to the instrument: with ``axis`` given they are taken about
    the head axis, not about the vertical.  The edge points may then carry a
    third coordinate; without it they are taken at the height of the scanner.
    """
    verts = voxel_vertices(lo, h)
    o = np.asarray(origin, dtype=float)
    ok = np.ones(verts.shape[0], dtype=bool)
    w = verts - o
    if is_vertical(axis):
        with np.errstate(divide="ignore", invalid="ignore"):
            zen = np.arccos(w[:, :, 2] / np.sqrt((w ** 2).sum(axis=2)))
        ok &= np.all((zen >= zenith_min) & (zen <= zenith_max), axis=1)
        if left_edge is not None and right_edge is not None:
            ok &= np.all(in_sector(verts[:, :, 0], verts[:, :, 1], o, left_edge, right_edge), axis=1)
        return ok
    pol = axis_angles(w, axis)[0]
    ok &= np.all((pol >= zenith_min) & (pol <= zenith_max), axis=1)
    if left_edge is not None and right_edge is not None:
        ok &= np.all(in_sector_about(w, o, left_edge, right_edge, axis), axis=1)
    return ok


def covered_matlab(centre, h, origin, left_edge, right_edge, pattern, zenith_min):
    """Coverage test of the MATLAB code (check_covered.m and the max_zen test).

    Signed distances d1, d2 of the voxel centre (x, y) from the two edge
    lines; the voxel is dropped when they have the same sign or either is
    within half a voxel of its line.  The test cannot tell the scanned sector
    from the one mirrored behind the scanner.  If an edge point equals the
    scanner position (the documented way to declare a 360-degree scan) the
    distances are NaN and every voxel passes.
    """
    c = np.asarray(centre, dtype=float)
    o = np.asarray(origin, dtype=float)
    ok = np.ones(c.shape[0], dtype=bool)
    if left_edge is not None and right_edge is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            x1, y1 = left_edge
            x2, y2 = right_edge
            n1 = np.sqrt((x1 - o[0]) ** 2 + (y1 - o[1]) ** 2)
            d1 = ((c[:, 0] - o[0]) * (y1 - o[1]) - (c[:, 1] - o[1]) * (x1 - o[0])) / n1
            n2 = np.sqrt((x2 - o[0]) ** 2 + (y2 - o[1]) ** 2)
            d2 = ((c[:, 0] - o[0]) * (y2 - o[1]) - (c[:, 1] - o[1]) * (x2 - o[0])) / n2
            bad = ((d1 > 0) & (d2 > 0)) | ((d1 < 0) & (d2 < 0)) | (np.abs(d1) < h / 2) | (np.abs(d2) < h / 2)
        ok &= ~bad
    if pattern == "spherical":
        ok &= ~(view_zenith(c, o) < zenith_min)
    return ok
