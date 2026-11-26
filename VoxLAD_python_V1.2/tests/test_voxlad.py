import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import voxlad
from voxlad import geometry as geo, lad as ladm, leafangle as la, traverse, legacy
from voxlad.model import Scan, Settings, compute
import simulate

HERE = os.path.dirname(__file__)
EXAMPLE = os.path.join(HERE, "..", "examples", "shea_tree")
have_example = os.path.isfile(os.path.join(EXAMPLE, "a1s1_t10.asc"))


# ---------------------------------------------------------------- geometry
def test_projected_area_far_voxel_matches_parallel_projection():
    # far away the cross-section of a cube is h^2 (|ux| + |uy| + |uz|)
    rng = np.random.default_rng(0)
    lo = rng.uniform(-1, 1, (200, 3)) * 500 + np.array([3000.0, 2000.0, 1000.0])
    h = 0.4
    a, nf = geo.projected_area(lo, h, [0, 0, 0])
    c = lo + h / 2
    u = c / np.linalg.norm(c, axis=1)[:, None]
    assert np.allclose(a, h * h * np.abs(u).sum(axis=1), rtol=1e-3)
    assert np.all(nf == 3)


def test_projected_area_face_rules():
    lo = np.array([[-0.1, 30.0, -0.05]])              # scanner between two pairs of faces
    exact, nf = geo.projected_area(lo, 0.3, [0, 0, 0], "exact")
    matlab, _ = geo.projected_area(lo, 0.3, [0, 0, 0], "matlab")
    assert nf[0] == 1
    assert np.isclose(exact[0], 0.09 * (30.15 / 30.0) ** 2, rtol=1e-4)
    assert matlab[0] > exact[0]
    lo = np.array([[5.0, 30.0, 4.0]])                 # ordinary corner view: both rules agree
    assert np.isclose(geo.projected_area(lo, 0.3, [0, 0, 0], "exact")[0][0],
                      geo.projected_area(lo, 0.3, [0, 0, 0], "matlab")[0][0], rtol=1e-12)


def test_matlab_round():
    assert np.array_equal(geo.matlab_round([0.5, 1.5, 2.5, 2.4999, 3.0]), [1, 2, 3, 2, 3])


def test_sector_and_coverage():
    o = [0.0, 0.0, 0.0]
    # sector from west (-1, 0) clockwise through north to east (1, 0)
    assert geo.in_sector(0.0, 5.0, o, [-1, 0], [1, 0])
    assert not geo.in_sector(0.0, -5.0, o, [-1, 0], [1, 0])
    lo = np.array([[-0.15, 10.0, 0.0], [-0.15, -10.3, 0.0], [10.0, -0.15, 0.0]])
    cov = geo.covered(lo, 0.3, o, [-1, 0], [1, 0])
    assert cov.tolist() == [True, False, False]       # third voxel straddles the east edge
    # the MATLAB test cannot tell the scanned sector from its mirror image behind the scanner
    centre = np.array([[0.0, 10.0, 0.0], [0.0, -10.0, 0.0]])
    assert geo.covered_matlab(centre, 0.3, o, [-1, 1], [1, 1], "uniform", 0.0).tolist() == [True, True]
    assert geo.covered(centre - 0.15, 0.3, o, [-1, 1], [1, 1]).tolist() == [True, False]
    # zenith limits: voxel straight overhead is outside a 30-130 degree field of view
    lo = np.array([[-0.15, -0.15, 10.0], [10.0, 0.0, 0.0]])
    assert geo.covered(lo, 0.3, o, zenith_min=np.radians(30), zenith_max=np.radians(130)).tolist() == [False, True]


def test_pulse_count_patterns():
    c = np.array([[0.0, 10.0, 10.0]])                 # zenith 45 degrees
    a = np.array([1.0])
    st = 1e-3
    n_uni = geo.pulse_count(a, c, [0, 0, 0], "uniform", st, st)[0]
    n_sph = geo.pulse_count(a, c, [0, 0, 0], "spherical", st, st)[0]
    assert np.isclose(n_uni, 1.0 / (200 * 1e-6), rtol=1e-3)
    assert np.isclose(n_sph / n_uni, np.sqrt(2), rtol=1e-3)   # 1/sin(45 deg) more pulses


# ---------------------------------------------------------------- leaf angle, LAD
def test_g_function_spherical_is_half_and_symmetric():
    th = np.radians(np.linspace(1, 179, 60))
    G = la.g_function(th, la.leaf_angle_table("spherical"))
    assert np.all(np.abs(G - 0.5) < 0.003)
    assert np.allclose(G, G[::-1])
    for name, g in la.LEAF_ANGLE_TABLES.items():
        assert abs(sum(g) - 1) < 1e-6, name


def test_g_matlab_formula_goes_complex_below_the_horizon():
    g = la.leaf_angle_table("erectophile")
    lit = la.g_function(np.radians([80.0, 94.0, 96.0, 120.0]), g, folded=False)
    assert np.all(lit.imag[:2] == 0) and np.all(lit.imag[2:] != 0)
    fold = la.g_function(np.radians([80.0, 94.0, 96.0, 120.0]), g)
    assert np.allclose(lit.real, fold, rtol=1e-12)    # the real part is right; only the imaginary part is spurious


def test_g_horizontal_view_is_continuous():
    g = la.leaf_angle_table("erectophile")
    a, b, c = la.g_function(np.array([np.pi / 2 - 1e-6, np.pi / 2, np.pi / 2 + 1e-6]), g)
    assert abs(a - b) < 1e-5 and abs(b - c) < 1e-5


def test_jensen_scaling_with_voxel_size():
    # the paper's term scales with h^2; the MATLAB constant does not
    for h in (0.1, 0.3, 1.0):
        dbar = h / 1.3
        assert np.isclose(ladm.jensen_delta(h, dbar, "paper") / h ** 2, ladm.jensen_delta(1.0, 1 / 1.3, "paper"))
    assert np.isclose(ladm.jensen_delta(0.1, 0.1 / 1.3, "matlab"), ladm.jensen_delta(1.0, 1 / 1.3, "matlab"))
    assert np.isclose(ladm.jensen_delta(0.3, 0.2, "matlab"), ladm.jensen_delta(0.3, 0.2, "matlab-scaled"))
    assert abs(ladm.jensen_delta(0.3, 0.2, "paper") / ladm.jensen_delta(0.3, 0.2, "matlab") - 1) < 0.05


def test_invert_lad_cases():
    dbar, G = 0.25, 0.5
    delta = ladm.jensen_delta(0.3, dbar, "paper")
    lad, st = ladm.invert_lad(dbar, delta, G, np.array([900.0, 1000.0, 0.0, 10.0]), np.array([100.0, 0.0, 50.0, 990.0]))
    assert st.tolist() == [ladm.OK, ladm.OK, ladm.NO_GAP, ladm.NO_SOLUTION]
    assert lad[1] == 0 and np.isnan(lad[2]) and np.isnan(lad[3])
    k = lad[0] * G                                     # the solution satisfies ln T = -k d + k^2 delta
    assert np.isclose(np.log(0.9), -k * dbar + k * k * delta)
    lad0, _ = ladm.invert_lad(dbar, 0.0, G, 900.0, 100.0)
    assert np.isclose(lad0, -np.log(0.9) / (G * dbar))  # no correction: Beer-Lambert with the mean path
    assert ladm.invert_lad(dbar, delta, G, 10.0, 990.0, "zero")[0] == 0
    assert np.isclose(ladm.invert_lad(dbar, delta, G, 10.0, 990.0, "beer")[0], -np.log(0.01) / (G * dbar))


def test_invert_lad_matlab_special_cases():
    dbar, G = 0.25, 0.5
    d = ladm.jensen_delta(0.3, dbar, "matlab")
    f = lambda ng, ni, g=G: float(ladm.invert_lad_matlab(dbar, d, g, ng, ni))
    assert np.isclose(f(900, 100), (dbar - np.sqrt(dbar ** 2 + 4 * d * np.log(0.9))) / (2 * d * G))
    assert f(0, 50) == 0            # log(0): complex -> 0
    assert f(-5, 50) == 0           # log of a negative number -> 0
    assert f(10, 990) == 0          # negative discriminant -> 0
    assert np.isnan(f(0, 0))        # 0/0
    assert f(-50, 50) == -np.inf    # Ng < 0 and Ng + Ni = 0
    assert f(-50, 10) < 0           # Ng < 0 and Ng + Ni < 0: a negative LAD is printed
    assert f(900, 100, 0.5 + 0.01j) == 0


def test_exact_path_length_inversion():
    h = 0.3
    u = np.array([[0.0, 1.0, 0.0], [0.5, 1.0, 0.3], [1.0, 1.0, 1.0]])
    u /= np.linalg.norm(u, axis=1)[:, None]
    # face-on: every ray crosses h, so T = exp(-k h) exactly
    assert np.isclose(ladm.transmission_exact(2.0, u[0], h), np.exp(-2.0 * h))
    # mean path = volume / cross-section: dT/dk at 0 is -h/(a+b+c)
    eps = 1e-6
    for ui in u:
        slope = (ladm.transmission_exact(eps, ui, h) - 1.0) / eps
        assert np.isclose(-slope, h / np.abs(ui).sum(), rtol=1e-4)
    # round trip, including optical depths where Eq. 8 has no solution
    for k in (0.5, 5.0, 40.0):
        T = ladm.transmission_exact(k, u, h)
        lad, st = ladm.invert_lad_exact(u, h, 0.5, T * 1e6, (1 - T) * 1e6)
        assert np.allclose(lad * 0.5, k, rtol=1e-6) and np.all(st == ladm.OK)
    # the published fit agrees with the exact solution within 2% up to K = k*h = 1.5
    for ui in u:
        dbar = h / np.abs(ui).sum()
        for K in (0.3, 1.0, 1.5):
            T = float(ladm.transmission_exact(K / h, ui, h))
            lad, st = ladm.invert_lad(dbar, ladm.jensen_delta(h, dbar, "paper"), 1.0, T, 1 - T)
            assert st == ladm.OK and abs(lad * h / K - 1) < 0.02
    lad, st = ladm.invert_lad_exact(u, h, 0.5, np.array([0.0, 5.0, 5.0]), np.array([3.0, 0.0, 5.0]))
    assert st.tolist() == [ladm.NO_GAP, ladm.OK, ladm.OK] and lad[1] == 0 and lad[2] > 0


# ---------------------------------------------------------------- engines
def _random_scene(n=30000, seed=2):
    rng = np.random.default_rng(seed)
    xyz = np.column_stack([rng.uniform(-3, 3, n), rng.uniform(8, 14, n), rng.uniform(-1, 4, n)])
    cls = rng.integers(1, 4, n)
    return xyz, cls


def _grid(shape=(20, 20, 17), origin=(-3.00003, 8.00003, -1.00003), h=0.3):
    nvox = int(np.prod(shape))
    return np.array(origin), h, shape, np.arange(nvox, dtype=np.int32), nvox


def test_traversal_equals_brute_force():
    xyz, cls = _random_scene()
    origin, h, shape, lookup, nvox = _grid()
    o = np.array([0.2, -1.0, 0.5])
    nb = traverse.count_blocked(xyz, o, origin, h, shape, lookup, nvox)
    rng = np.random.default_rng(5)
    for flat in rng.choice(nvox, 60, replace=False):
        ijk = np.array(np.unravel_index(flat, shape))
        l = origin + ijk * h
        d = xyz - o
        ta, tb = (l - o) / d, (l + h - o) / d
        t0 = np.minimum(ta, tb).max(axis=1)
        t1 = np.maximum(ta, tb).min(axis=1)
        inside = np.all((xyz >= l) & (xyz < l + h), axis=1)
        assert nb[flat] == np.sum((t0 < t1) & (t0 >= 1) & ~inside)


def test_traversal_numpy_fallback_identical():
    xyz, cls = _random_scene(8000)
    origin, h, shape, lookup, nvox = _grid()
    for o in ([0.2, -1.0, 0.5], [0.0, 11.0, 1.5], [-10.0, 11.0, 30.0]):   # outside, inside and above the grid
        a = traverse.count_blocked(xyz, o, origin, h, shape, lookup, nvox, use_numba=False)
        if traverse.HAVE_NUMBA:
            b = traverse.count_blocked(xyz, o, origin, h, shape, lookup, nvox)
            assert np.array_equal(a, b)
        assert a.sum() > 0


def test_returns_outside_the_grid_still_block():
    origin, h, shape, lookup, nvox = np.array([0.0, 10.0, 0.0]), 1.0, (1, 1, 1), np.zeros(1, np.int32), 1
    o = np.array([0.5, 0.0, 0.5])
    pts = np.array([[0.5, 5.0, 0.5],      # in front of the voxel, on the line of sight: blocks
                    [0.5, 10.5, 0.5],     # inside: a hit, not a blocked pulse
                    [0.5, 15.0, 0.5],     # behind
                    [3.0, 5.0, 0.5]])     # off to the side
    assert traverse.count_blocked(pts, o, origin, h, shape, lookup, nvox)[0] == 1
    assert traverse.count_hits(pts, np.array([1, 1, 1, 1]), origin, h, shape, lookup, nvox)[0].tolist() == [1, 0, 0]


@pytest.mark.skipif(not legacy.HAVE_NUMBA, reason="needs numba")
def test_matlab_engine_agrees_with_traversal_in_corner_views():
    # scanner outside all three slabs of every voxel: both engines must agree
    # (grid shifted so that no return lies on a voxel face)
    xyz, cls = _random_scene(20000)
    origin, h, shape, lookup, nvox = _grid()
    o = np.array([-20.0, -15.0, -9.0])
    nb = traverse.count_blocked(xyz, o, origin, h, shape, lookup, nvox)
    ni = traverse.count_hits(xyz, cls, origin, h, shape, lookup, nvox)
    ijk = np.stack(np.unravel_index(np.arange(nvox), shape), axis=1)
    sel = np.random.default_rng(1).choice(nvox, 400, replace=False)
    m = legacy.count_matlab(xyz, cls, origin + ijk[sel] * h, h, o)
    assert np.array_equal(m[:, 0], nb[sel])
    assert np.array_equal(m[:, 1:], ni[sel])


# ---------------------------------------------------------------- end to end
@pytest.mark.parametrize("pattern,box", [("spherical", ([5.0, -1.0, 5.0], [7.0, 1.0, 7.0])),
                                         ("uniform", ([5.0, -1.0, 5.0], [7.0, 1.0, 7.0])),
                                         ("spherical", ([9.0, -1.0, -1.0], [11.0, 1.0, 1.0]))])
def test_recovers_known_lad_from_simulated_scan(pattern, box):
    rng = np.random.default_rng(11)
    lo, hi = box
    true_lad = 1.0
    g = la.leaf_angle_table("spherical")
    step = np.radians(0.04)
    scans = []
    for pos in ([0.0, 0.0, 0.0], [0.0, 12.0, 0.5]):
        xyz, _ = simulate.scan_box(pos, lo, hi, true_lad, g, step, step, rng, pattern)
        scans.append(Scan(xyz, np.ones(len(xyz), int), pos, pattern, step))
    res = compute(scans, Settings(voxel_size=0.5, bounds=lo + hi, voxels="all", leaf_angle="spherical"))
    exact = compute(scans, Settings(voxel_size=0.5, bounds=lo + hi, voxels="all", leaf_angle="spherical",
                                    jensen_fit="exact"))
    assert abs(np.nanmean(exact.lad) - true_lad) < 0.03
    assert np.nanmax(np.abs(exact.lad - res.lad)) < 0.05      # LAD*G*h = 0.25 here: the two agree closely
    assert res.lad.size == 64
    assert np.all(res.flag <= voxlad.FLAG_LOW_SAMPLING)
    good = res.flag == voxlad.FLAG_OK
    assert good.sum() >= 60
    # each voxel is estimated from a few thousand pulses; the mean over the box must be close
    assert abs(np.mean(res.lad[good]) - true_lad) < 0.03
    assert np.std(res.lad[good]) < 0.15


def test_low_sampling_and_unresolved_flags():
    # a dense slab in front of a second box: the far box is poorly sampled
    rng = np.random.default_rng(3)
    g = la.leaf_angle_table("spherical")
    step = np.radians(0.06)
    pos = [0.0, 0.0, 0.5]
    front, _ = simulate.scan_box(pos, [-1, 5, 0], [1, 6, 1], 12.0, g, step, step, rng)
    scan = Scan(front, np.ones(len(front), int), pos, "spherical", step)
    res = compute([scan], Settings(voxel_size=0.5, bounds=[-1, 5, 0, 1, 7, 1], voxels="all"))
    k = res.ijk[:, 1]
    assert np.all(res.flag[k >= 2] >= voxlad.FLAG_LOW_SAMPLING)     # behind the slab
    assert np.all(res.explored[k == 0] > 90)                  # a few pulses cross a dense neighbour first
    assert np.nanmax(res.lad[k == 0]) > 5


def test_unresolved_voxels_and_fallbacks():
    # a very dense box seen obliquely from one side: front voxels transmit a few percent only
    rng = np.random.default_rng(8)
    g = la.leaf_angle_table("spherical")
    step = np.radians(0.08)
    pos = [0.0, 0.0, 4.0]
    lo, hi = [2.0, 4.0, 0.0], [3.0, 5.0, 1.0]
    xyz, _ = simulate.scan_box(pos, lo, hi, 14.0, g, step, step, rng)
    scan = Scan(xyz, np.ones(len(xyz), int), pos, "spherical", step)
    kw = dict(voxel_size=0.5, bounds=lo + hi, voxels="all")
    res = compute([scan], Settings(**kw))
    unres = res.flag == voxlad.FLAG_UNRESOLVED
    assert unres.any() and np.all(np.isnan(res.lad[unres]))
    assert not np.any(np.isnan(res.lad[res.flag <= voxlad.FLAG_LOW_SAMPLING]))
    zero = compute([scan], Settings(unresolved="zero", **kw))
    assert np.all(zero.lad[unres] == 0) and np.array_equal(zero.flag, res.flag)
    beer = compute([scan], Settings(unresolved="beer", **kw))
    st = res.per_scan["status"][:, 0]
    assert np.all(np.isfinite(beer.lad[st == ladm.NO_SOLUTION])) and np.all(np.isnan(beer.lad[st == ladm.NO_GAP]))
    exact = compute([scan], Settings(jensen_fit="exact", **kw))
    assert np.all(np.isfinite(exact.lad[st == ladm.NO_SOLUTION]))
    assert np.all(exact.lad[st == ladm.NO_SOLUTION] > beer.lad[st == ladm.NO_SOLUTION])   # Beer-Lambert is a lower bound
    # a second scan from the other side resolves voxels the first one cannot
    pos2 = [5.0, 9.0, 4.0]
    xyz2, _ = simulate.scan_box(pos2, lo, hi, 14.0, g, step, step, rng)
    both = compute([scan, Scan(xyz2, np.ones(len(xyz2), int), pos2, "spherical", step)], Settings(**kw))
    assert (both.flag == voxlad.FLAG_UNRESOLVED).sum() < unres.sum()
    picked = both.best_scan[np.isfinite(both.lad)] - 1
    rows = np.flatnonzero(np.isfinite(both.lad))
    assert np.all(both.per_scan["status"][rows, picked] == ladm.OK)


# ---------------------------------------------------------------- files, config, CLI
def test_config_template_is_valid(tmp_path):
    p = tmp_path / "t.toml"
    p.write_text(voxlad.TEMPLATE)
    settings, specs, out = voxlad.load_config(str(p))
    assert settings.voxel_size == 0.3 and len(specs) == 2 and specs[0]["pattern"] == "spherical"


def test_config_errors(tmp_path):
    p = tmp_path / "bad.toml"
    p.write_text('[grid]\nvoxel_size = 0.3\nvoxle = 1\n[[scan]]\nfile="a"\n')
    with pytest.raises(voxlad.InputError, match="unknown setting"):
        voxlad.load_config(str(p))
    p.write_text('[grid]\nvoxel_size = 0.3\n')
    with pytest.raises(voxlad.InputError, match="no \\[\\[scan\\]\\]"):
        voxlad.load_config(str(p))


def test_text_and_asc_readers_and_class_mapping(tmp_path):
    rng = np.random.default_rng(0)
    g = la.leaf_angle_table("spherical")
    step = np.radians(0.1)
    xyz, _ = simulate.scan_box([0, 0, 0], [4, -1, -1], [6, 1, 1], 1.0, g, step, step, rng)
    asc = tmp_path / "s.asc"
    with open(asc, "w") as fh:
        fh.write("0 0 0\n0 0\n0 0\nd 0.1\n%d\n" % len(xyz))
        np.savetxt(fh, np.column_stack([xyz, np.ones(len(xyz))]), fmt="%.4f %.4f %.4f %d")
    csv = tmp_path / "s.csv"
    np.savetxt(csv, np.column_stack([xyz, np.full(len(xyz), 7), np.ones(len(xyz))]), fmt="%.4f,%.4f,%.4f,%d,%d",
               header="x,y,z,class,rn", comments="")
    cfg = tmp_path / "c.toml"
    cfg.write_text('''
[grid]
voxel_size = 0.5
bounds = [4, -1, -1, 6, 1, 1]
[model]
leaf_angle = "spherical"
[classes]
leaf = [1, 7]
[scanner]
pattern = "spherical"
[[scan]]
file = "s.asc"
[[scan]]
file = "s.csv"
position = [0, 0, 0]
angular_step_deg = 0.1
columns = ["x", "y", "z", "class", "return_number"]
skip_rows = 1
[output]
prefix = "out/run"
npz = true
matlab_text = true
''')
    res = voxlad.run(str(cfg))
    # the same returns read from both formats: identical counts for the two scans
    assert np.array_equal(res.per_scan["ni_leaf"][:, 0], res.per_scan["ni_leaf"][:, 1])
    assert np.array_equal(res.per_scan["nb"][:, 0], res.per_scan["nb"][:, 1])
    for suffix in ("_voxels.csv", "_scans.csv", "_summary.txt", "_grid.npz", "_matlab.txt"):
        assert (tmp_path / "out" / ("run" + suffix)).is_file()
    tab = np.genfromtxt(tmp_path / "out" / "run_voxels.csv", delimiter=",", names=True)
    assert np.allclose(tab["lad"], res.lad, rtol=1e-5, equal_nan=True)
    z = np.load(tmp_path / "out" / "run_grid.npz")
    assert z["lad"].shape == res.shape
    # unknown class code
    cfg.write_text(cfg.read_text().replace("leaf = [1, 7]", "leaf = [1]"))
    with pytest.raises(voxlad.InputError, match="class codes \\[7\\]"):
        voxlad.run(str(cfg))


def test_cli(tmp_path, capsys):
    from voxlad.cli import main
    assert main(["template", str(tmp_path / "x.toml")]) == 0
    assert main(["run", str(tmp_path / "missing.toml")]) == 2
    assert "not found" in capsys.readouterr().err


# ---------------------------------------------------------------- instrument not level
def _angle(a, b):
    return np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b))))))


def test_pulse_count_follows_the_scanner_axis():
    """Turning scanner, voxels and axis together changes nothing; the vertical formulas are the special case."""
    rng = np.random.default_rng(8)
    h, o = 0.4, np.array([0.5, -1.0, 1.2])
    lo = rng.uniform(-15, 15, (400, 3))
    c = lo + h / 2
    ap, _ = geo.projected_area(lo, h, o)
    st = np.radians(0.05)
    n0 = geo.pulse_count(ap, c, o, "spherical", st, st)
    assert np.array_equal(n0, geo.pulse_count(ap, c, o, "spherical", st, st, axis=[0, 0, 2.0]), equal_nan=True)
    R = simulate.tilt_matrix(17.0, 250.0)
    a = R[:, 2]
    v = c - o
    d = np.linalg.norm(v, axis=1)
    pol, _ = geo.axis_angles(v, a)
    n1 = geo.pulse_count(ap, c, o, "spherical", st, st, axis=a)
    expect = ap / (2 * d * np.tan(st / 2) * 2 * d * np.sin(pol) * np.tan(st / 2))
    assert np.allclose(n1, geo.matlab_round(expect))
    assert np.allclose(geo.axis_distance(v, a), d * np.sin(pol))
    # the error made by taking the axis as vertical is the ratio of the sines
    zen = np.arccos(v[:, 2] / d)
    far = np.sin(zen) > 0.3
    assert np.allclose((n0 / n1)[far], (np.sin(pol) / np.sin(zen))[far], rtol=0.01)
    assert voxlad.orientation.tilt_effect(1.0, 60.0) == pytest.approx(0.0304, abs=1e-3)
    # limits of the scan are taken about the axis: a voxel 40 degrees from a lying axis is inside 30..130
    lying = np.array([1.0, 0.0, 0.0])
    lo1 = np.array([[10 * np.cos(np.radians(40)), 0.0, 10 * np.sin(np.radians(40))]]) - 0.2
    assert geo.covered(lo1, 0.4, np.zeros(3), None, None, np.radians(30), np.radians(130), axis=lying)[0]
    assert not geo.covered(lo1, 0.4, np.zeros(3), None, None, np.radians(55), np.radians(130), axis=lying)[0]
    assert geo.covered(lo1, 0.4, np.zeros(3), None, None, np.radians(30), np.radians(130))[0]      # zenith 50 degrees


@pytest.mark.parametrize("kw,tilt,position_error", [
    ({}, 2.3, (0.0, 0.0, 0.0)),                                                  # planes through the axis
    ({"cone_deg": 1.0, "origin_offset": 0.004}, 0.4, (0.0, 0.0, 0.0)),           # lines as a polygon mirror draws them
    ({"head_range": (20.0, 90.0)}, 7.0, (0.0, 0.0, 0.0)),                        # a sector scan
    ({}, 88.0, (0.0, 0.0, 0.0)),                                                 # scanner on a tilt mount
    ({}, 1.2, (0.15, -0.13, 0.0)),                                               # declared position 20 cm off
])
def test_axis_found_from_the_scan_lines(kw, tilt, position_error):
    from voxlad.orientation import estimate_axis
    rng = np.random.default_rng(4)
    step = np.radians(0.25)
    pos = np.array([1.0, 2.0, 0.3])
    xyz, axis = simulate.scan_room(pos, simulate.tilt_matrix(tilt, 130.0, 40.0), step, rng, **kw)
    est = estimate_axis(np.round(xyz, 3), pos + np.array(position_error), step)        # a file written in millimetres
    assert est.ok, est.reason
    assert _angle(est.axis, axis) < 0.01
    assert est.tilt_deg == pytest.approx(tilt, abs=0.01) and est.bearing_deg == pytest.approx(130.0, abs=1.5)
    assert est.se_deg < 0.01 and est.rms_steps < 0.1
    assert est.cone_deg == pytest.approx(kw.get("cone_deg", 0.0), abs=0.02)
    # the lines meet at the scanner: an error of the declared position across the axis is found
    assert np.allclose(est.offset[:2], -np.array(position_error[:2]), atol=0.01)
    assert est.step_between == pytest.approx(step, rel=0.01) and est.step_along == pytest.approx(step, rel=0.01)


def test_axis_is_left_alone_when_the_file_cannot_tell():
    from voxlad.orientation import estimate_axis, resolve_axis
    rng = np.random.default_rng(5)
    step = np.radians(0.25)
    pos = np.array([1.0, 2.0, 0.3])
    xyz, axis = simulate.scan_room(pos, simulate.tilt_matrix(3.0, 20.0), step, rng)
    est = estimate_axis(xyz[rng.permutation(len(xyz))], pos, step)             # the order of acquisition is lost
    assert not est.ok and "scan lines" in est.reason
    assert not estimate_axis(xyz[:500], pos, step).ok
    cls = np.ones(len(xyz), int)
    shuffled = Scan(xyz[rng.permutation(len(xyz))], cls, pos, "spherical", step)
    info = resolve_axis(shuffled)
    assert info["axis"] is None and info["source"] == "vertical" and "could not be found" in info["text"]
    # ... unless a registration matrix is at hand
    M = np.eye(4)
    M[:3, :3] = simulate.tilt_matrix(3.0, 20.0)
    M[:3, 3] = pos
    info = resolve_axis(Scan(shuffled.xyz, cls, pos, "spherical", step, matrix=M))
    assert info["source"] == "matrix" and _angle(info["axis"], axis) < 1e-6 and info["tilt_deg"] == pytest.approx(3.0)
    # settings: a vector, "vertical", MATLAB mode; a wrong declared step and a wrong given axis are reported
    good = Scan(xyz, cls, pos, "spherical", step)
    assert resolve_axis(good)["source"] == "estimated" and resolve_axis(good)["warnings"] == []
    assert resolve_axis(good, matlab=True)["axis"] is None
    assert resolve_axis(Scan(xyz, cls, pos, "spherical", step, axis="vertical"))["axis"] is None
    given = resolve_axis(Scan(xyz, cls, pos, "spherical", step, axis=[0, 0, 1]))
    assert given["axis"] is None and len(given["warnings"]) == 1 and "3.00 deg from the axis found" in given["warnings"][0]
    w = resolve_axis(Scan(xyz, cls, pos, "spherical", step * 1.05))["warnings"]
    assert len(w) == 1 and "angular step" in w[0] and "10% higher" in w[0]
    w = resolve_axis(Scan(xyz, cls, pos + [0.3, 0.0, 0.0], "spherical", step))["warnings"]
    assert len(w) == 1 and "0.30 m from the declared scanner position" in w[0]
    with pytest.raises(voxlad.InputError):
        Scan(xyz, cls, pos, "spherical", step, axis="sideways")
    with pytest.raises(voxlad.InputError):
        Scan(xyz, cls, pos, "spherical", step, axis="matrix")
    with pytest.raises(voxlad.InputError):
        Scan(xyz, cls, pos, "spherical", step, matrix=np.ones((4, 4)))


def test_tilted_scanner_gives_the_same_lad_as_a_level_one():
    """A box of leaves of known density seen steeply by a scanner leaning 5 degrees towards it."""
    rng = np.random.default_rng(12)
    g = la.leaf_angle_table("spherical")
    st = np.radians(0.04)
    lo, hi = [3.0, -1.0, 6.0], [5.0, 1.0, 8.0]
    settings = Settings(voxel_size=0.5, bounds=lo + hi, voxels="all", leaf_angle="spherical", jensen_fit="exact")
    R = simulate.tilt_matrix(5.0, 90.0)
    xyz, _ = simulate.scan_box([0, 0, 0], lo, hi, 1.0, g, st, st, rng, rotation=R)
    cls = np.ones(len(xyz), int)
    res = compute([Scan(xyz, cls, [0, 0, 0], "spherical", st, name="leaning")], settings)
    info = res.scan_axes[0]
    assert info["source"] == "estimated" and info["tilt_deg"] == pytest.approx(5.0, abs=0.02)
    assert info["bearing_deg"] == pytest.approx(90.0, abs=0.5)
    assert abs(np.nanmean(res.lad) - 1.0) < 0.05 and np.all(res.flag == voxlad.FLAG_OK)
    assert "head axis: tilted 5.00 deg" in voxlad.summary_text(res)
    # the polar angle about the axis differs from the zenith angle by up to the tilt; G keeps the zenith angle
    dz = res.per_scan["polar_deg"][:, 0] - res.per_scan["zenith_deg"][:, 0]
    assert 4.0 < np.abs(dz).max() <= 5.0 + 1e-6
    given = compute([Scan(xyz, cls, [0, 0, 0], "spherical", st, axis=R[:, 2])], settings)
    assert np.allclose(given.lad, res.lad, rtol=2e-3)
    # taken as vertical, the same returns give far too much leaf area on the side the scanner leans to
    wrong = compute([Scan(xyz, cls, [0, 0, 0], "spherical", st, axis="vertical")], settings)
    assert np.nanmean(wrong.lad) > 1.25
    # MATLAB mode keeps the vertical axis of the original program
    if legacy.HAVE_NUMBA:
        old = compute([Scan(xyz, cls, [0, 0, 0], "spherical", st)], Settings(voxel_size=0.5, compat="matlab"))
        assert old.scan_axes[0]["axis"] is None and "MATLAB mode" in old.scan_axes[0]["text"]


def test_axis_settings_in_the_configuration_file(tmp_path):
    rng = np.random.default_rng(6)
    step = 0.5
    pos = np.array([1.0, 2.0, 0.3])
    R = simulate.tilt_matrix(4.0, 300.0)
    xyz, axis = simulate.scan_room(pos, R, np.radians(step), rng)
    np.savetxt(tmp_path / "s.txt", np.column_stack([xyz, np.ones(len(xyz))]), fmt="%.3f %.3f %.3f %d")
    M = np.eye(4)
    M[:3, :3], M[:3, 3] = R, pos
    np.savetxt(tmp_path / "s_sop.dat", M, fmt="%.9f")
    base = ('[grid]\nvoxel_size = 1.0\nbounds = [3, 3, 2, 6, 6, 5]\nvoxels = "all"\n[scanner]\nformat = "text"\n'
            'pattern = "spherical"\nangular_step_deg = %g\n%s[[scan]]\nfile = "s.txt"\n%s[output]\nprefix = "o/t"\n')

    def load(scanner="", scan="position = [1.0, 2.0, 0.3]\n"):
        (tmp_path / "t.toml").write_text(base % (step, scanner, scan))
        _, specs, _ = voxlad.load_config(str(tmp_path / "t.toml"))
        return voxlad.load_scan(specs[0])
    sc = load()
    assert sc.axis == "auto" and sc.matrix is None
    assert load('axis = "vertical"\n').axis == "vertical"
    assert np.allclose(load('axis = [0.0, 0.3, 3.0]\n').axis, np.array([0, 0.3, 3.0]) / np.hypot(0.3, 3.0))
    m = load(scan='matrix_file = "s_sop.dat"\n')                       # the position comes with the matrix
    assert np.allclose(m.position, pos) and np.allclose(m.matrix, M, atol=1e-8)
    inline = ", ".join("%.9f" % v for v in M.ravel())
    assert np.allclose(load(scan="position = [1.0, 2.0, 0.3]\nmatrix = [%s]\naxis = \"matrix\"\n" % inline).matrix, M, atol=1e-8)
    for bad in ('axis = "up"\n', 'axis = [1, 2]\n', 'axis = "matrix"\n'):
        with pytest.raises(voxlad.InputError):
            load(bad)
    with pytest.raises(voxlad.InputError):
        load(scan='position = [1.0, 2.0, 0.3]\nmatrix = [1, 2, 3]\n')
    with pytest.raises(voxlad.InputError):
        load(scan='position = [1.0, 2.0, 0.3]\nmatrix_file = "nothing.dat"\n')
    # a run: the axis is found and reported, the scan table carries the polar angle
    sc = load()
    res = voxlad.run(str(tmp_path / "t.toml"))
    assert res.scan_axes[0]["source"] == "estimated" and _angle(res.scan_axes[0]["axis"], axis) < 0.02
    assert "polar_deg" in open(tmp_path / "o" / "t_scans.csv").readline()
    assert "head axis: tilted 4.00 deg" in open(tmp_path / "o" / "t_summary.txt").read()
    from voxlad import scaninfo
    text = scaninfo.describe(sc)
    assert "head axis      tilted 4.00 deg" in text and "spacing of the returns: 0.5 deg between lines" in text
    # a scanner lying on its side needs the limits of the scan as points in space
    lying = Scan(sc.xyz, sc.cls, pos, "spherical", sc.step_az, None, (0.0, 5.0), (5.0, 0.0), axis=[1.0, 0.0, 0.05])
    with pytest.raises(voxlad.InputError, match="x, y, z"):
        compute([lying], Settings(voxel_size=1.0, bounds=[3, 3, 2, 6, 6, 5], voxels="all"))
    lying = Scan(sc.xyz, sc.cls, pos, "spherical", sc.step_az, None, (0.0, 5.0, 3.0), (5.0, 0.0, -3.0), axis=[1.0, 0.0, 0.05])
    compute([lying], Settings(voxel_size=1.0, bounds=[3, 3, 2, 6, 6, 5], voxels="all"))


# ---------------------------------------------------------------- 3-D viewer
def _viewer_payload(path):
    import base64, json, re
    page = open(path, encoding="utf-8").read()
    m = re.search(r'<script type="application/json" id="voxlad-data">(.*?)</script>', page, re.S)
    d = json.loads(m.group(1))
    arr = {k: np.frombuffer(base64.b64decode(v["b64"]), dtype=np.dtype(v["dtype"]).newbyteorder("<"))
           for k, v in d["arrays"].items()}
    return page, d, arr


def _small_result():
    rng = np.random.default_rng(5)
    lo, hi = [5.0, -1.0, 5.0], [7.0, 1.0, 7.0]
    step = np.radians(0.1)
    scans = []
    for n, pos in enumerate(([0.0, 0.0, 0.0], [0.0, 12.0, 0.5])):
        xyz, _ = simulate.scan_box(pos, lo, hi, 1.0, la.leaf_angle_table("spherical"), step, step, rng, "spherical")
        scans.append(Scan(xyz, np.ones(len(xyz), int), pos, "spherical", step, name="s%d <b>" % (n + 1)))
    return compute(scans, Settings(voxel_size=0.5, bounds=lo + hi, voxels="all", leaf_angle="spherical"))


def test_viewer_holds_the_result(tmp_path):
    res = _small_result()
    res.lad[3] = np.nan                                   # a voxel without a value keeps its place
    path = voxlad.write_viewer(res, str(tmp_path / "out" / "box_viewer.html"), csv="box_voxels.csv")
    page, d, arr = _viewer_payload(path)
    assert page.lstrip().startswith("<!DOCTYPE html>")
    assert d["scans"][0]["name"] == "s1 <b>" and "s1 <b>" not in page     # names cannot inject markup into the page
    assert "http://" not in page and "https://" not in page          # nothing is fetched from the network
    assert d["title"] == "box" and d["n"] == res.lad.size == 64 and d["voxel_size"] == 0.5
    assert d["csv"] == "box_voxels.csv" and d["version"] == voxlad.__version__
    assert [s["name"] for s in d["scans"]] == res.scan_names
    assert np.allclose([s["position"] for s in d["scans"]], [[0, 0, 0], [0, 12, 0.5]])
    assert np.array_equal(arr["ijk"].reshape(-1, 3), res.ijk)
    assert np.array_equal(arr["lad"], res.lad.astype(np.float32), equal_nan=True) and np.isnan(arr["lad"][3])
    assert np.array_equal(arr["flag"], res.flag) and np.array_equal(arr["scan"], res.best_scan)
    rows = np.arange(64)
    assert np.array_equal(arr["nt"], res.per_scan["nt"][rows, res.best_scan - 1])
    # voxel centres as the viewer works them out
    centre = np.asarray(d["origin"]) + (arr["ijk"].reshape(-1, 3) + 0.5) * d["voxel_size"]
    assert np.allclose(centre, res.lo + 0.25)


def test_viewer_from_voxel_table(tmp_path):
    res = _small_result()
    files = voxlad.write_results(res, str(tmp_path / "box"))
    assert files[-1].endswith("box_viewer.html")
    _, d0, a0 = _viewer_payload(files[-1])
    os.remove(files[-1])
    path = voxlad.write_viewer_from_csv(str(tmp_path / "box_voxels.csv"))
    assert path == files[-1]
    _, d1, a1 = _viewer_payload(path)
    assert d1["voxel_size"] == d0["voxel_size"] and np.allclose(d1["origin"], d0["origin"], atol=1e-4)
    assert d1["scans"] == []                                   # the voxel table does not hold the scanner positions
    for k in ("ijk", "flag", "scan", "nt", "nb", "nil", "niw", "nin"):
        assert np.array_equal(a0[k], a1[k]), k
    assert np.allclose(a0["lad"], a1["lad"], atol=1e-6, equal_nan=True)
    assert voxlad.write_results(res, str(tmp_path / "plain"), viewer=False)[-1].endswith("_summary.txt")
    (tmp_path / "other.csv").write_text("a,b\n1,2\n")
    with pytest.raises(voxlad.InputError):
        voxlad.write_viewer_from_csv(str(tmp_path / "other.csv"))


def test_viewer_opens_only_when_asked(tmp_path, monkeypatch, capsys):
    import shutil
    import webbrowser
    from voxlad import viewer
    from voxlad.cli import main
    opened = []
    monkeypatch.setattr(webbrowser, "open", lambda url, *a, **k: opened.append(url) or True)
    monkeypatch.setattr(viewer, "can_open_browser", lambda: True)
    res = _small_result()
    np.savetxt(tmp_path / "a.txt", np.column_stack([res.lo[:40] + 0.2, np.ones(40)]), fmt="%.3f")
    cfg = tmp_path / "t.toml"
    text = ('[grid]\nvoxel_size = 0.5\n[scanner]\nformat = "text"\npattern = "spherical"\nangular_step_deg = 0.5\n'
            '[[scan]]\nfile = "a.txt"\nposition = [0, 0, 0]\n[output]\nprefix = "out/t"\n')
    cfg.write_text(text)
    assert main(["run", str(cfg), "-q", "--no-view"]) == 0 and opened == []
    html = tmp_path / "out" / "t_viewer.html"
    assert html.is_file()
    assert main(["run", str(cfg), "-q"]) == 0                              # the default: the viewer opens
    assert len(opened) == 1 and opened[0].startswith("file:") and opened[0].endswith("t_viewer.html")
    cfg.write_text(text + "open_viewer = false\n")
    assert main(["run", str(cfg), "-q"]) == 0 and len(opened) == 1
    cfg.write_text(text + "viewer = false\n")
    os.remove(html)
    assert main(["run", str(cfg), "-q"]) == 0 and len(opened) == 1 and not html.exists()
    # voxlad view: from the configuration file, a voxel table or the page itself; the page is written
    # from the table when it is missing, older than the table, or with --rebuild
    table = tmp_path / "out" / "t_voxels.csv"
    assert main(["view", str(cfg), "--no-open"]) == 0 and html.is_file() and len(opened) == 1
    assert "wrote" in capsys.readouterr().out
    assert main(["view", str(cfg)]) == 0 and len(opened) == 2
    assert main(["view", str(table)]) == 0 and len(opened) == 3
    assert main(["view", str(html)]) == 0 and len(opened) == 4
    assert "wrote" not in capsys.readouterr().out                           # an up-to-date page is left alone
    assert main(["view", str(html), "--rebuild", "--no-open"]) == 0 and "wrote" in capsys.readouterr().out
    t = os.path.getmtime(html)
    os.utime(table, (t + 5, t + 5))                                         # the table was changed after the page
    assert main(["view", str(table), "--no-open"]) == 0 and "wrote" in capsys.readouterr().out
    assert main(["view", str(tmp_path / "nothing.html")]) == 2
    (tmp_path / "alone.html").write_text("<html></html>")
    assert main(["view", str(tmp_path / "alone.html")]) == 0 and len(opened) == 5      # any page opens
    assert main(["view", str(tmp_path / "alone.html"), "--rebuild"]) == 2              # but has no table to rebuild from
    capsys.readouterr()
    opened.pop()
    # from Python: run() opens nothing unless asked; view() writes a page and opens it
    cfg.write_text(text)
    voxlad.run(str(cfg))
    assert len(opened) == 4
    voxlad.run(str(cfg), view=True)
    assert len(opened) == 5 and opened[-1].endswith("t_viewer.html")
    voxlad.run(str(cfg), write=False, view=True)
    assert len(opened) == 6
    p = voxlad.view(res, str(tmp_path / "v.html"), open=False)
    assert os.path.isfile(p) and len(opened) == 6
    p = voxlad.view(res)
    assert os.path.isfile(p) and len(opened) == 7
    os.remove(p)
    # no graphical session: nothing is started and the run says where the page is
    monkeypatch.setattr(viewer, "can_open_browser", lambda: False)
    assert main(["run", str(cfg)]) == 0 and len(opened) == 7
    assert "open " in capsys.readouterr().out
    shutil.rmtree(tmp_path / "out")


# ---------------------------------------------------------------- the worked example
@pytest.mark.skipif(not have_example, reason="example scans not present")
def test_shea_tree_example(tmp_path):
    res = voxlad.run(os.path.join(EXAMPLE, "shea_tree.toml"), write=False)
    assert res.lad.size == 4402
    assert abs(res.total_leaf_area() - 37.595) < 0.01
    assert np.all(res.flag <= voxlad.FLAG_LOW_SAMPLING)


REFERENCE = os.path.join(HERE, "..", "validation", "octave", "reference", "thin100_dbg.txt")


@pytest.mark.skipif(not (have_example and legacy.HAVE_NUMBA and os.path.isfile(REFERENCE)),
                    reason="example scans, Octave reference or numba not present")
def test_matlab_mode_pipeline_against_octave():
    """Whole pipeline in MATLAB mode (grid, occupancy map, scan limits, counts)
    against voxlad_v7.m run in GNU Octave on the scans thinned 1:100.

    validation/compare_octave.py shows that all counts are identical when the
    voxel corners are built exactly as Octave builds its loop variable.  Here
    the corners are origin + k * voxel_size; a few returns lying exactly on a
    voxel face then fall on the other side, which moves single counts.
    """
    import math
    settings, specs, _ = voxlad.load_config(os.path.join(EXAMPLE, "shea_tree.toml"))
    k = 100
    thin = [Scan(s.xyz[::k], s.cls[::k], s.position, "spherical", float("%.8f" % (s.step_az * math.sqrt(k))),
                 None, s.left_edge, s.right_edge, np.radians(32.0), np.pi, s.name)       # TLS_type = 2, max_zen = 32
            for s in voxlad.load_scans(specs)]
    res = compute(thin, Settings(voxel_size=0.3, z_range=(-1, 12), leaf_angle="erectophile", compat="matlab"))
    ref = np.loadtxt(REFERENCE)
    key = {tuple(np.round(r, 3)): i for i, r in enumerate(res.lo)}
    assert len(key) == len(ref) // 2 == 2016                       # the same voxels are evaluated
    same, worst = 0, 0
    for r in ref:
        i, s = key[tuple(np.round(r[:3], 3))], int(r[3]) - 1
        mine = np.array([res.per_scan[c][i, s] for c in ("nt", "nb", "ni_leaf", "ni_wood", "ni_noise")])
        same += np.array_equal(mine, r[8:13])
        worst = max(worst, int(np.abs(mine - r[8:13]).max()))
    assert same >= 0.99 * len(ref) and worst <= 1
