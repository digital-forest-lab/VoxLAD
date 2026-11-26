"""Leaf angle distributions and the leaf projection function G(theta, theta_L).

A leaf angle distribution is given as the fraction of leaf area in each of
``n`` equal-width leaf inclination classes between 0 and 90 degrees (the
original MATLAB code uses nine 10-degree classes centred on 5, 15, ... 85
degrees).  The projection function follows Ross (1981), as in the MATLAB code.
"""
from __future__ import annotations

import numpy as np

# Nine-class tables copied digit for digit from voxlad_v7.m (lines 59-63).
LEAF_ANGLE_TABLES = {
    "erectophile": [0.001813, 0.014469, 0.038908, 0.072257, 0.11037,
                    0.148837, 0.182804, 0.208301, 0.222241],
    "planophile": [0.217993772, 0.205566574, 0.181747777, 0.149126381,
                   0.111758486, 0.07430429, 0.040992495, 0.016051798,
                   0.002458428],
    "spherical": [0.013816774, 0.044002464, 0.072956086, 0.099621579,
                  0.1232069, 0.143184018, 0.158760891, 0.169409487,
                  0.175041802],
    "plagiophile": [0.007084426, 0.053135816, 0.126252094, 0.19247555,
                    0.220395847, 0.197361602, 0.133842925, 0.059766886,
                    0.009684853],
    "uniform": [0.111111111] * 9,
}
_ALIASES = {"erecto": "erectophile", "plano": "planophile",
            "spheri": "spherical", "plagio": "plagiophile"}

# Beam-size correction factors H per leaf inclination class, copied from
# voxlad_v7.m (line 106).  They were derived for the Optech ILRIS-3D at about
# 20 m range and one leaf size; they do not apply to other instruments.
H_ILRIS3D_20M = [0.445537, 0.598981, 0.669567, 0.710171, 0.733597,
                 0.747112, 0.75491, 0.759409, 0.762004]


def leaf_angle_table(spec):
    """Return class fractions as a float array from a name or a sequence."""
    if isinstance(spec, str):
        key = _ALIASES.get(spec.lower(), spec.lower())
        if key not in LEAF_ANGLE_TABLES:
            raise ValueError(
                f"unknown leaf angle distribution '{spec}'; choose one of "
                f"{sorted(LEAF_ANGLE_TABLES)} or give a list of class fractions")
        return np.asarray(LEAF_ANGLE_TABLES[key], dtype=float)
    g = np.asarray(spec, dtype=float).ravel()
    if g.size < 1 or np.any(g < 0) or not np.all(np.isfinite(g)):
        raise ValueError("leaf angle class fractions must be finite and >= 0")
    return g


def class_centres(n):
    """Mid inclination (radians) of n equal classes spanning 0-90 degrees."""
    # degrees first, then *(pi/180): the same arithmetic as deg2rad(5), deg2rad(15)...
    return ((np.arange(n) + 0.5) * (90.0 / n)) * (np.pi / 180.0)


def g_function(theta, g, folded=True):
    """Leaf projection function G for view zenith angle(s) ``theta`` (rad).

    folded=True (default) uses G(theta) = G(pi - theta), so directions below
    the horizontal are handled.  folded=False evaluates the MATLAB expression
    literally and returns a complex array: for theta > 90 deg + theta_L the
    arccos argument is below -1, which is what makes the MATLAB code write
    LAD = 0 for voxels well below scanner height.
    """
    theta = np.asarray(theta, dtype=float)
    g = np.asarray(g, dtype=float)
    tl = class_centres(g.size)
    if folded:
        th = np.where(theta > np.pi / 2, np.pi - theta, theta)
        G = np.zeros(th.shape)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            for gi, tli in zip(g, tl):
                cc = np.cos(th) * np.cos(tli)
                ss = np.sin(th) * np.sin(tli)
                # cot(th)cot(tl) = cc/ss ; >= 1 exactly when th <= pi/2 - tl
                c = np.clip(np.where(ss > 0, cc / np.where(ss > 0, ss, 1.0), np.inf), -1.0, 1.0)
                x = np.arccos(c)
                # algebraically identical to cc*(1 + 2*(tan(x) - x)/pi), but
                # stable at theta = 90 degrees
                s_else = cc * (1.0 - 2.0 * x / np.pi) + (2.0 / np.pi) * ss * np.sin(x)
                G += gi * np.where(th <= np.pi / 2 - tli, cc, s_else)
        return G
    # literal MATLAB evaluation (complex where MATLAB goes complex)
    G = np.zeros(theta.shape, dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cot_t = 1.0 / np.tan(theta)
        for gi, tli in zip(g, tl):
            c = cot_t * (1.0 / np.tan(tli))
            first = theta <= np.pi / 2 - tli
            s1 = np.cos(theta) * np.cos(tli)
            inrange = np.abs(c) <= 1.0
            xr = np.arccos(np.where(inrange, c, 0.0))
            s = np.where(first, s1, s1 * (1.0 + (2.0 * (np.tan(xr) - xr) / np.pi))).astype(complex)
            cplx = (~first) & (~inrange) & np.isfinite(c)
            if np.any(cplx):
                xc = np.arccos(c[cplx].astype(complex))
                s[cplx] = s1[cplx] * (1.0 + (2.0 * (np.tan(xc) - xc) / np.pi))
            if gi != 0.0:
                G = G + gi * s
    return G
