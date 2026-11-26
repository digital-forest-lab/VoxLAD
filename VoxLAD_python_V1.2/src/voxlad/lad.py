"""Leaf area density from pulse counts (Beland et al. 2014, Eqs. 1 and 7-9)."""
from __future__ import annotations

import numpy as np

JENSEN_FITS = ("paper", "matlab", "matlab-scaled", "none", "exact")


def jensen_delta(voxel_size, d_bar, fit="paper"):
    """Path-length variability term, in m^2 (written delta_d in the MATLAB code).

    With K = LAD*G, the model is  ln T = -K*d_bar + K^2*delta.  In the paper's
    notation delta = Delta_d * Delta_h, where Eq. 7 gives the ratio
    Delta_d/Delta_h, so delta scales with the square of the voxel size.

    fit = "paper"          Eq. 7: 0.07162*(1 - exp(-9.536*(Dh/d_bar - 1))) * Dh^2
    fit = "matlab"         0.006239*(1 - exp(-10.38*(Dh/d_bar - 1))); the MATLAB
                           constant, in m^2, right only for 0.3 m voxels
    fit = "matlab-scaled"  the MATLAB fit multiplied by (Dh/0.3)^2
    fit = "none"           0: Beer-Lambert with the mean path length (Eq. 2)

    ("exact" is not a value of delta: see invert_lad_exact.)
    """
    h = float(voxel_size)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r = h / np.asarray(d_bar, dtype=float) - 1.0
        if fit == "paper":
            return 0.07162 * h * h * (1.0 - np.exp(-9.536 * r))
        if fit == "matlab":
            return 0.006239 * (1.0 - np.exp(-10.38 * r))
        if fit == "matlab-scaled":
            return 0.006239 * (h / 0.3) ** 2 * (1.0 - np.exp(-10.38 * r))
        if fit == "none":
            return np.zeros_like(r)
    raise ValueError(f"jensen_fit must be one of {JENSEN_FITS}, not '{fit}'")


# status of one voxel seen by one scan
OK = 0
NOT_COVERED = 1      # outside the scan's azimuth sector or zenith range
TOO_CLOSE = 2        # scanner inside or too near the voxel, or voxel at the scan pole
NO_PULSES = 3        # theoretical pulse count rounds to zero
OCCLUDED = 4         # no pulse reached the voxel (Nt - Nb <= 0)
NO_GAP = 5           # every pulse that entered was intercepted (Ng <= 0)
NO_SOLUTION = 6      # transmission too low: Eq. 8 has no real solution
SCAN_STATUS = {OK: "ok", NOT_COVERED: "not_covered", TOO_CLOSE: "too_close",
               NO_PULSES: "no_pulses", OCCLUDED: "occluded", NO_GAP: "no_gap",
               NO_SOLUTION: "no_solution"}


def invert_lad(d_bar, delta, G, ng, ni_leaf_star, unresolved="nan"):
    """Solve Eq. 8 for LAD.  All arguments broadcast.

    Returns (lad, status) where status is OK, NO_GAP or NO_SOLUTION.  The
    caller is expected to have excluded uncovered and occluded cases.

    unresolved = "nan"   leave NO_GAP / NO_SOLUTION voxels undefined
               = "zero"  write 0 for them (what the MATLAB code does)
               = "beer"  for NO_SOLUTION use -ln(T)/(G*d_bar), a lower bound;
                         NO_GAP stays undefined
    """
    d_bar = np.asarray(d_bar, dtype=float)
    delta = np.asarray(delta, dtype=float)
    G = np.asarray(G, dtype=float)
    ng = np.asarray(ng, dtype=float)
    nil = np.asarray(ni_leaf_star, dtype=float)
    shape = np.broadcast(d_bar, delta, G, ng, nil).shape
    lad = np.full(shape, np.nan)
    status = np.full(shape, OK, dtype=np.int8)
    with np.errstate(divide="ignore", invalid="ignore"):
        no_gap = ng <= 0
        T = np.where(no_gap, np.nan, ng / (ng + nil))
        lnT = np.log(T)
        disc = d_bar ** 2 + 4.0 * delta * lnT
        bad = (disc < 0) & ~no_gap
        root = np.sqrt(np.where(disc >= 0, disc, np.nan))
        # (d - sqrt(d^2 + 4*delta*lnT)) / (2*delta*G), written so that it stays
        # accurate as delta -> 0 (it then tends to -lnT/(G*d))
        val = -2.0 * lnT / (G * (d_bar + root))
        lad = np.where(no_gap | bad, np.nan, val)
        beer = -lnT / (G * d_bar)
    status = np.where(no_gap, NO_GAP, np.where(bad, NO_SOLUTION, OK)).astype(np.int8)
    if unresolved == "zero":
        lad = np.where(no_gap | bad, 0.0, lad)
    elif unresolved == "beer":
        lad = np.where(bad, beer, lad)
    elif unresolved != "nan":
        raise ValueError("unresolved must be 'nan', 'zero' or 'beer'")
    return lad, status


def invert_lad_matlab(d_bar, delta, G, ng, ni_leaf_star):
    """LAD exactly as the MATLAB script prints it, including its special cases.

    MATLAB evaluates
        LAD = (d_bar - sqrt(d_bar^2 + 4*delta_d*log(Ng/(Ng+Ni))))/(2*delta_d*G)
    in complex arithmetic and then sets LAD = 0 whenever the result is not
    real.  ``G`` may be complex (see leafangle.g_function with folded=False).
    """
    d_bar = np.asarray(d_bar, dtype=float)
    delta = np.asarray(delta, dtype=float)
    Gc = np.asarray(G, dtype=complex)
    ng = np.asarray(ng, dtype=float)
    nil = np.asarray(ni_leaf_star, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r = ng / (ng + nil)
        g_complex = Gc.imag != 0
        Gr = Gc.real
        lnr = np.log(np.where(r > 0, r, 1.0))
        disc = d_bar ** 2 + 4.0 * delta * lnr
        num = d_bar - np.sqrt(np.where(disc >= 0, disc, 0.0))
        real_val = num / (2.0 * delta * Gr)
        # complex result (printed as 0): log of a negative number or of zero
        # (sqrt(-Inf)), negative discriminant, or complex G
        goes_complex = (np.isfinite(r) & (r < 0)) | (r == 0) | (np.isfinite(r) & (r > 0) & (disc < 0))
        lad = np.where(np.isnan(r), np.where(g_complex, 0.0, np.nan),
                       np.where(goes_complex | g_complex, 0.0, real_val))
        # Ng < 0 with Ng + Ni = 0: log(-Inf) = Inf + pi*i, and the square root of
        # (Inf + b*i) is Inf with zero imaginary part, so the result is real -Inf
        lad = np.where(np.isneginf(r) & ~g_complex, -np.inf / np.sign(delta * Gr), lad)
    return lad


# --------------------------------------------------------------------------
# Exact treatment of the path lengths for parallel rays (an extension; not in
# the 2014 paper).
#
# For parallel rays of direction u = (a, b, c) (absolute values of a unit
# vector) through a cube of side h, the rays whose path inside the cube is at
# least l are those that cross both the cube and the cube shifted by l*u.  The
# intersection is a box of sides h - l*a, h - l*b, h - l*c, whose cross-section
# along u is
#     A(l) = a(h-lb)(h-lc) + b(h-la)(h-lc) + c(h-la)(h-lb)
#          = h^2*s1 - 2*h*s2*l + 3*s3*l^2,   s1=a+b+c, s2=ab+ac+bc, s3=abc
# for 0 <= l <= L = h/max(a, b, c).  Hence path lengths have a density that is
# linear in l, (2*h*s2 - 6*s3*l)/A(0), plus a share A(L)/A(0) of rays with the
# maximal length L (those going through two opposite faces), and the mean
# transmission  T(k) = E[exp(-k*l)]  has a closed form.  It is solved for
# k = LAD*G without any fitted coefficient, and has a solution for every
# transmission between 0 and 1.
# --------------------------------------------------------------------------
def _chord_terms(u_abs, h):
    u = np.abs(np.asarray(u_abs, dtype=float))
    a, b, c = u[..., 0], u[..., 1], u[..., 2]
    s1 = a + b + c
    s2 = a * b + a * c + b * c
    s3 = a * b * c
    A0 = h * h * s1
    L = h / np.maximum(np.maximum(a, b), c)
    p0 = 2.0 * h * s2 / A0
    p1 = -6.0 * s3 / A0
    m = np.clip((A0 - 2.0 * h * s2 * L + 3.0 * s3 * L * L) / A0, 0.0, 1.0)
    return p0, p1, m, L, h / s1


def transmission_exact(k, u_abs, voxel_size):
    """Mean transmission of parallel rays of direction u through a cube filled
    with a turbid medium of extinction coefficient k = LAD*G (per metre)."""
    p0, p1, m, L, _ = _chord_terms(u_abs, voxel_size)
    return _t_exact(np.asarray(k, dtype=float), p0, p1, m, L)


def _t_exact(k, p0, p1, m, L):
    x = k * L
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        e = np.exp(-x)
        small = x < 1e-2
        xs = np.where(small, x, 1.0)                    # avoid 0/0 in the unused branch
        xl = np.where(small, 1.0, x)
        phi0 = np.where(small, 1.0 - xs / 2 + xs ** 2 / 6 - xs ** 3 / 24, -np.expm1(-xl) / xl)
        phi1 = np.where(small, 0.5 - xs / 3 + xs ** 2 / 8 - xs ** 3 / 30 + xs ** 4 / 144,
                        (1.0 - np.exp(-xl) * (1.0 + xl)) / xl ** 2)
    return p0 * L * phi0 + p1 * L * L * phi1 + m * e


def invert_lad_exact(u_abs, voxel_size, G, ng, ni_leaf_star, d_bar=None):
    """LAD from the exact path-length distribution (see above).

    u_abs : (n, 3) direction scanner -> voxel centre (unit vector; signs ignored)
    d_bar : mean path length of the model, h^3/A_proj.  The pulses of a scanner
            diverge, which shortens the mean path slightly compared with
            parallel rays (by about h/d for a voxel seen face-on at range d).
            When d_bar is given, the path lengths are scaled so that their
            mean equals it.
    Returns (lad, status) with status OK or NO_GAP.
    """
    G = np.asarray(G, dtype=float)
    ng = np.asarray(ng, dtype=float)
    nil = np.asarray(ni_leaf_star, dtype=float)
    p0, p1, m, L, dbar = _chord_terms(u_abs, voxel_size)
    no_gap = ng <= 0
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(no_gap, 0.5, ng / (ng + nil))      # placeholder where undefined
        lo = -np.log(T) / dbar                          # Beer-Lambert with the mean path: a lower bound
        hi = 2.0 * lo + 1e-9
        for _ in range(80):                             # grow the bracket
            need = _t_exact(hi, p0, p1, m, L) > T
            if not np.any(need):
                break
            hi = np.where(need, hi * 2.0, hi)
        for _ in range(70):                             # bisection: T(k) decreases with k
            mid = 0.5 * (lo + hi)
            above = _t_exact(mid, p0, p1, m, L) > T
            lo = np.where(above, mid, lo)
            hi = np.where(above, hi, mid)
        k = 0.5 * (lo + hi)
        if d_bar is not None:
            k = k * (dbar / np.asarray(d_bar, dtype=float))
        lad = np.where(no_gap, np.nan, np.where(T >= 1.0, 0.0, k / G))
    return lad, np.where(no_gap, NO_GAP, OK).astype(np.int8)
