"""Path-length (Jensen) term: Monte Carlo, closed form, published fit, MATLAB constant.

Parallel rays through a cube of side h; chord lengths d_i by the slab method.
Truth: T(k) = mean(exp(-k d_i)) with k = LAD*G.  Each method then recovers k
from T; the tables give recovered / true.
"""
import numpy as np

from _common import voxlad  # noqa: F401
from voxlad.lad import invert_lad, invert_lad_exact, jensen_delta, transmission_exact

rng = np.random.default_rng(1)


def chords(u, h, n=600000):
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    a = np.cross(u, [0, 0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0, 0])
    a /= np.linalg.norm(a)
    b = np.cross(u, a)
    R = h * np.sqrt(3) / 2
    p = rng.uniform(-R, R, (n, 1)) * a + rng.uniform(-R, R, (n, 1)) * b + h / 2 - 10 * h * u
    with np.errstate(divide="ignore", invalid="ignore"):
        ta, tb = (0 - p) / u, (h - p) / u
    t0, t1 = np.minimum(ta, tb).max(1), np.maximum(ta, tb).min(1)
    return (t1 - t0)[t1 > t0], u


views = {"face-on": [0, 1, 0.0], "3 deg": [0.05, 1, 0.02], "10 deg": [0.18, 1, 0.03],
         "oblique": [0.5, 1, 0.3], "diagonal": [1, 1, 1]}
print("A. Recovered / true LAD; K = LAD*G*h is the optical depth across one voxel side")
print(f"{'view':9s} {'h':>4s} {'K':>4s} {'T (MC)':>8s} {'T closed':>9s} {'exact':>7s} {'paper':>7s} {'matlab':>7s} {'none':>7s}")
for name, u in views.items():
    for h in (0.1, 0.3, 1.0):
        d, uu = chords(u, h)
        dbar = d.mean()
        for K in (0.5, 1.5):
            k = K / h
            T = float(np.mean(np.exp(-k * d)))
            Tc = float(transmission_exact(k, uu[None, :], h)[0])
            row = [float(invert_lad_exact(uu[None, :], h, 1.0, np.array([T]), np.array([1 - T]))[0][0]) / k]
            for fit in ("paper", "matlab", "none"):
                lad, st = invert_lad(dbar, jensen_delta(h, dbar, fit), 1.0, T, 1 - T)
                row.append(float(lad) / k if st == 0 else np.nan)
            print(f"{name:9s} {h:4.1f} {K:4.1f} {T:8.4f} {Tc:9.4f} " + " ".join(f"{v:7.3f}" for v in row))

print()
print("B. Published fit (Eq. 7 with Eq. 8) against the closed form, 4000 random directions, h = 0.3 m")
u = np.abs(rng.normal(size=(4000, 3)))
u /= np.linalg.norm(u, axis=1)[:, None]
h = 0.3
dbar = h / u.sum(1)
for K in (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
    T = transmission_exact(K / h, u, h)
    lad, st = invert_lad(dbar, jensen_delta(h, dbar, "paper"), 1.0, T, 1 - T)
    r = lad * h / K
    beer = -np.log(T) / dbar * h / K
    print(f"K={K:4.2f}  T {T.min():.3f}..{T.max():.3f}  solved {100 * np.mean(st == 0):5.1f}%  "
          f"recovered/true mean {np.nanmean(r):.3f} range {np.nanmin(r):.3f}..{np.nanmax(r):.3f}  "
          f"no correction {beer.min():.3f}..{beer.max():.3f}")

print()
print("C. Lowest transmission for which Eq. 8 has a real solution")
for r in (1.05, 1.1, 1.2, 1.3, 1.5, 1.732):
    a = 0.07162 * (1 - np.exp(-9.536 * (r - 1)))
    line = f"h/d_bar = {r:5.3f}: paper fit T > {np.exp(-(1 / r) ** 2 / (4 * a)):.3f}"
    for hh in (0.1, 0.3, 1.0):
        dm = 0.006239 * (1 - np.exp(-10.38 * (r - 1)))
        line += f"; MATLAB constant at {hh} m T > {np.exp(-(hh / r) ** 2 / (4 * dm)):.3f}"
    print(line)
