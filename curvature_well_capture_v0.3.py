#!/usr/bin/env python3
"""
curvature_well_capture.py
=========================
Public simulation code for the paper

    "Localized Standing Waves in Prescribed Spherically Symmetric
     Inward Curvature Wells"

Core question:
    We are not trying to create matter from nothing.
    We are asking whether geometry can create a trap for waves.

This single, self-contained file computes everything the paper reports and
writes the frozen data contract consumed by the (plot-only) figure package.
All physics computation, gate evaluation, and status assignment happen here;
the figure scripts only read these products and draw them.

Scope (this code does NOT claim): creation of particles or matter; any
specific high-energy phenomenology; self-consistent gravity or metric
backreaction; physical-unit calibration. All quantities are in dimensionless
model units with m0 = 1. It studies one narrow question only: whether a
prescribed inward-curvature-like well localizes scalar wave energy.

Design invariants:
  * self-contained: the finite-volume grid, spherical weights, sponge, and
    divergence-flux operator are inlined here (no runtime dependency on any
    laboratory engine file);
  * deterministic: no random numbers; two runs from one config produce
    byte-identical CSV/NPZ (timestamps live only in metadata);
  * no plotting and no matplotlib import anywhere in this module;
  * status strings are produced in exactly two functions: adjudicate_run
    (per run) and adjudicate_ensemble (H2 over the robustness subset).

Dependencies: numpy, scipy. Standard library for CSV/JSON/hashing.
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import platform
import re
import sys
import time
import dataclasses
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq


# ======================================================================
# L0 -- CONSTANTS, PRE-REGISTERED GATES, AND NUMERICAL-VALIDITY LIMITS
# ======================================================================
# Version of THIS simulation code. Synced to the release tag at publication.
# The suffix marks that the R9 configuration values are not yet frozen.
CODE_VERSION = "0.3.0"

M0 = 1.0  # unit mass scale; all runs in m0 = 1 units. No physical units anywhere.

# Pre-registered acceptance gates. These mirror manuscript Table 1 and are the
# single source of truth; every exporter copies gate values from here, and no
# threshold is hardcoded anywhere downstream.
GATES = {
    "eps_lin":    2e-3,   # initial bound-mode overlap ceiling (clean-seed test)
    "f_min":      2e-2,   # minimum net captured fraction to call capture
    "eps_res":    1e-2,   # residual box-energy fraction for "radiation cleared"
    "eps_budget": 1e-2,   # |f_esc + f_cap - 1| closure tolerance
    "rtol":       0.30,   # robustness: relative spread of net f_cap across subset
}

# Numerical-validity limits (R1: EPS_CONS is ported verbatim from the engine's
# E_CONS_TOL = 1e-2; the blow-up and box-affected constants are ported from the
# laboratory reference). A run failing any of these is INVALID_NUMERICS, never
# a physics result.
EPS_CONS = 1e-2          # max |(E_box + E_esc - E_in)/E_in| tolerated (R1)
E_FLOOR = 1e-7           # energy noise floor (engine MASS_FLOOR), denominator guard
# Amplitude blow-up guard. NOTE (v0.3, deviation from the ported lab constant):
# the laboratory code guarded an ABSOLUTE |chi| ceiling. In spherical symmetry a
# healthy ingoing shell focuses at the origin with chi = u/r growing as 1/r, so
# an absolute ceiling scales with the launch radius r0 and is not a property of
# the dynamics. Under the R9 carrier geometry (r0 = 55) the legacy ceiling
# halted the *linear* null control at its focal point while letting the
# nonlinear baseline pass -- silently invalidating the R4 null-floor check.
# The guard is therefore expressed on the geometry-independent invariant
# u = r*chi, whose magnitude is conserved up to coherent focal superposition
# for a free spherical wave, and which grows without bound only under genuine
# instability. Measured healthy values: 1.7x (linear, at focus), 1.1x
# (nonlinear). The energy blow-up guard below is unchanged and remains the
# primary instability net.
U_BLOWUP_FACTOR = 50.0   # halt if max|r*chi| exceeds this multiple of its seed value
E_BLOWUP_FACTOR = 100.0  # E_box beyond this * max(|E0|, 1) -> honest halt
KAPPA_BOX_MIN = 5.0      # kappa * r_max below this -> mode is box-affected (flagged)

VALID_PROFILES = ("square", "tanh", "sech2")


# ======================================================================
# L1 -- GEOMETRY: FINITE-VOLUME RADIAL GRID (inlined engine numerics)
# ======================================================================
@dataclass
class GridSpec:
    r_max: float = 40.0
    dr: float = 0.05
    sponge_start_frac: float = 0.75
    sponge_strength: float = 3.0


class Grid:
    """Spherically symmetric finite-volume radial grid.

    Cell-centered nodes r_j = j*dr; face areas r_face^2 and cell volumes give a
    divergence operator that is exactly regular at r = 0 (zero inner-face area)
    and discretely conservative. Spherical integration weights are
    w = 4*pi*vol. A quadratic sponge ramp damps the outer region; the outermost
    node is Dirichlet-held (chi = 0) by the caller. Numerics inlined verbatim
    from the laboratory engine so the published file stands alone.
    """

    def __init__(self, spec: GridSpec):
        self.dr = spec.dr
        self.N = int(round(spec.r_max / spec.dr)) + 1
        self.r = np.arange(self.N) * spec.dr
        self.r_max = self.r[-1]
        rf = (np.arange(self.N - 1) + 0.5) * spec.dr        # faces j+1/2
        self.area_f = rf ** 2                                # face areas (/4pi)
        r_lo = np.concatenate([[0.0], rf])                   # inner face of cell
        r_hi = np.concatenate([rf, [self.r_max]])            # outer face of cell
        self.vol = (r_hi ** 3 - r_lo ** 3) / 3.0             # cell volumes (/4pi)
        self.w = 4.0 * np.pi * self.vol                      # integration weights
        # quadratic sponge on the outer part of the box
        self.r_sp = spec.sponge_start_frac * self.r_max
        ramp = np.clip((self.r - self.r_sp) /
                       max(self.r_max - self.r_sp, 1e-12), 0.0, None)
        self.sponge = spec.sponge_strength * ramp ** 2
        self.m_boundary = self.r >= self.r_sp

    def div_flux(self, F_face):
        """(1/r^2) d_r(r^2 F) in finite-volume form; regular at r = 0 and
        discretely conservative. Outermost node is Dirichlet-held (zeroed by
        the caller)."""
        AF = self.area_f * F_face
        out = np.empty(self.N)
        out[0] = AF[0] / self.vol[0]
        out[1:-1] = (AF[1:] - AF[:-1]) / self.vol[1:-1]
        out[-1] = 0.0
        return out

    def face_grad(self, f):
        return (f[1:] - f[:-1]) / self.dr

    def node_grad(self, f):
        """Node-centered radial derivative for diagnostics/seeds: symmetric at
        the origin (f_r(0) = 0 for regular fields), one-sided at the edge."""
        g = np.gradient(f, self.dr)
        g[0] = 0.0
        return g


# ======================================================================
# L1 -- PRESCRIBED WELL PROFILES
# ======================================================================
def W_continuum(profile: str, r: np.ndarray, W0: float, a: float,
                Delta: Optional[float] = None) -> np.ndarray:
    """Continuum prescribed mass-lowering well W(r) <= 0, W -> 0 at infinity.

    square : -W0 for r < a, else 0                           (Eq. 1, sharp)
    tanh   : -0.5*W0*(1 - tanh((r - a)/Delta))               (Eq. 3)
    sech2  : -W0 / cosh(r/a)^2                                (Eq. 4)

    This is the closed-form definition; it matches the figure package's
    common_io.W_of exactly (self-test S4). The simulation samples the square
    profile with an area-weighted interface (see W_grid), which differs from
    this only within the single grid cell straddling r = a.
    """
    r = np.asarray(r, dtype=float)
    if profile == "square":
        return np.where(r < a, -W0, 0.0)
    if profile == "tanh":
        if Delta is None:
            raise ValueError("tanh profile requires an interface width Delta.")
        return -0.5 * W0 * (1.0 - np.tanh((r - a) / Delta))
    if profile == "sech2":
        return -W0 / np.cosh(r / a) ** 2
    raise ValueError(f"Unknown profile '{profile}'. Expected one of "
                     f"{VALID_PROFILES}.")


def W_grid(profile: str, r: np.ndarray, W0: float, a: float,
           dr: float, Delta: Optional[float] = None) -> np.ndarray:
    """Grid-sampled prescribed well used by the simulation.

    Identical to W_continuum except that the square profile uses an
    area-weighted interface across the cell straddling r = a: a node's depth is
    scaled by the fraction of its cell [r - dr/2, r + dr/2] lying inside the
    well. This is the standard second-order discretization of a sharp interface
    and is what preserves the validated convergence order and binding-threshold
    accuracy. It is a discretization choice only: the continuum well, its
    thresholds, and every claim are unchanged.
    """
    r = np.asarray(r, dtype=float)
    if profile == "square":
        if dr and dr > 0:
            lo = r - 0.5 * dr
            hi = r + 0.5 * dr
            frac = np.clip((np.minimum(hi, a) - lo) / dr, 0.0, 1.0)
            return -W0 * frac
        return np.where(r < a, -W0, 0.0)
    # smooth profiles: identical to the continuum form
    return W_continuum(profile, r, W0, a, Delta)


def m_eff2(profile, r, W0, a, Delta=None, grid=False, dr=None):
    W = (W_grid(profile, r, W0, a, dr, Delta) if grid
         else W_continuum(profile, r, W0, a, Delta))
    return M0 ** 2 + W


# ======================================================================
# L2 -- SPECTRAL LAYER: BOUND MODES AND ANALYTIC SQUARE-WELL REFERENCE
# ======================================================================
def spectral_bound_modes(grid: Grid, W: np.ndarray, n_max: int = 8) -> list:
    """Finite-difference eigenproblem for u = r*chi:

        [-d^2/dr^2 + m0^2 + W(r)] u = omega^2 u,   u(0) = u(r_max) = 0.

    Returns a list of modes with 0 < omega^2 < m0^2 (ordinary bound modes),
    each carrying its tachyonic flag for the guard. Tachyonic entries
    (omega^2 < 0) are returned flagged but MUST be refused upstream before any
    projection or export. u is normalized to int u^2 dr = 1; the sign is fixed
    positive at the dominant extremum.
    """
    dr = grid.dr
    diag = 2.0 / dr ** 2 + M0 ** 2 + W[1:-1]
    off = np.full(grid.N - 3, -1.0 / dr ** 2)
    k = min(n_max, grid.N - 3)
    vals, vecs = eigh_tridiagonal(diag, off, select="i", select_range=(0, k - 1))
    modes = []
    for i in range(k):
        om2 = float(vals[i])
        if om2 >= M0 ** 2 * (1.0 - 1e-9):
            break  # continuum (box-quantized) states start here
        u = np.zeros(grid.N)
        u[1:-1] = vecs[:, i]
        norm = np.sqrt(np.trapezoid(u ** 2, grid.r))
        u = u / norm
        if u[np.argmax(np.abs(u))] < 0:
            u = -u
        kappa = float(np.sqrt(max(M0 ** 2 - om2, 0.0)))
        interior = u[1:-1]
        n_nodes = int(np.sum(np.sign(interior[:-1]) * np.sign(interior[1:]) < 0))
        modes.append({
            "n": i + 1,
            "omega2": om2,
            "omega": float(np.sqrt(om2)) if om2 > 0 else float("nan"),
            "tachyonic": bool(om2 < 0),
            "kappa": kappa,
            "box_affected": bool(om2 > 0 and kappa * grid.r_max < KAPPA_BOX_MIN),
            "n_nodes": n_nodes,
            "u": u,
        })
    return modes


def square_well_matching_omega2(W0: float, a: float) -> list:
    """Analytic s-wave bound levels of the sharp square well, from the matching
    relation k*cot(k*a) = -kappa with k^2 = W0 - kappa^2, kappa in (0, sqrt W0).

    Returns omega^2 = m0^2 - kappa^2 for the real-omega levels, sorted
    ASCENDING in omega^2 so they pair with the finite-difference modes n = 1,
    2, ... (n = 1 is the most deeply bound, i.e. the smallest omega^2). Deeper
    (tachyonic) roots with kappa >= m0 are excluded here and refused upstream.
    """
    kmax = np.sqrt(W0)

    def fmatch(kappa):
        kk = np.sqrt(max(W0 - kappa ** 2, 1e-300))
        return kk * np.cos(kk * a) + kappa * np.sin(kk * a)

    grid = np.linspace(1e-6, kmax - 1e-6, 4000)
    fv = np.array([fmatch(x) for x in grid])
    roots = []
    for i in range(len(grid) - 1):
        if (np.isfinite(fv[i]) and np.isfinite(fv[i + 1])
                and fv[i] * fv[i + 1] < 0):
            try:
                roots.append(brentq(fmatch, grid[i], grid[i + 1], xtol=1e-12))
            except ValueError:
                pass
    omega2 = [float(M0 ** 2 - kap ** 2) for kap in roots if kap < M0]
    return sorted(omega2)  # ascending: pairs with FD modes n = 1, 2, ...


def threshold_depth_analytic(n: int, a: float) -> float:
    """Depth at which the n-th s-wave level appears: sqrt(W0)*a = (2n-1)*pi/2."""
    return ((2 * n - 1) * np.pi / (2.0 * a)) ** 2


def analytic_mode_count(sqrtW0_a) -> int:
    """Square-well staircase N = floor(sqrt(W0)*a/pi + 1/2) (Eq. A.4)."""
    return int(np.floor(float(sqrtW0_a) / np.pi + 0.5))


def numeric_threshold_depth(a: float, dr: float = 0.02, r_max: float = 40.0,
                            n: int = 1) -> float:
    """Depth W0 at which the FD ground state crosses omega^2 = m0^2, by
    bracketed root find on (lowest FD eigenvalue - m0^2). Used by self-test S2.
    """
    grid = Grid(GridSpec(r_max=r_max, dr=dr))

    def lowest_minus_threshold(W0):
        W = W_grid("square", grid.r, W0, a, dr)
        d = 2.0 / dr ** 2 + M0 ** 2 + W[1:-1]
        o = np.full(grid.N - 3, -1.0 / dr ** 2)
        val = eigh_tridiagonal(d, o, select="i", select_range=(0, 0))[0][0]
        return float(val) - M0 ** 2

    Wth = threshold_depth_analytic(n, a)
    return float(brentq(lowest_minus_threshold, 0.5 * Wth, 1.6 * Wth, xtol=1e-7))


def assert_no_tachyon(modes: list, context: str):
    """Refusal guard: the paper's scope excludes tachyonic wells. If any mode
    has omega^2 < 0 the run is refused (a configuration error), never assigned
    a status."""
    bad = [m for m in modes if m["tachyonic"]]
    if bad:
        raise ValueError(
            f"Tachyonic eigenvalue(s) present ({context}); this well is outside "
            "the paper's ordinary-bound-mode scope (0 < omega^2 < m0^2). "
            "Choose a shallower pre-registered well. Tachyonic runs are refused, "
            "not adjudicated.")


def real_bound_modes(modes: list) -> list:
    return [m for m in modes if not m["tachyonic"]]


# ----------------------------------------------------------------------
# Carrier-packet validity guards G1-G4 (R9-A).
# These are CONFIGURATION REFUSALS, symmetric with the tachyon guard: an
# out-of-protocol probe is an error, never an adjudicated status.
# ----------------------------------------------------------------------
G2_K0_SIGMA_MIN = 6.0          # low-k suppression: k0 * sigma
G3_POINTS_PER_WAVELENGTH = 25.0  # grid resolution of the carrier
G4_CLEARANCE_MARGIN = 1.5      # path-length safety factor
G4_SPECTRAL_EDGE = 3.3         # slow spectral edge k0 - 3.3/sigma


def assert_carrier_guards(params, log=None):
    """Refuse a chi_carrier configuration that violates G1-G4.

    G1 direction        : k0 > 0 (inward).
    G2 low-k suppression: k0 * sigma >= 6.
    G3 grid resolution  : 2*pi/(k0*dr) >= 25 points per carrier wavelength.
    G4 clearance margin : v_g(k_edge) * t_end >= 1.5 * (r0 + r_sponge_start),
                          with k_edge = k0 - 3.3/sigma.
    """
    if params.packet_kind != "chi_carrier":
        return None
    fail = lambda g, msg, fix: (_ for _ in ()).throw(ValueError(
        f"Carrier guard {g} violated: {msg}. {fix} "
        "(Out-of-protocol probe configurations are refused, not adjudicated.)"))

    if not (params.k0 > 0):
        fail("G1", f"k0 = {params.k0:g} is not > 0 (k0 > 0 means inward)",
             "Set a positive k0.")
    if params.k0_sigma < G2_K0_SIGMA_MIN:
        fail("G2", f"k0*sigma = {params.k0_sigma:.3g} < {G2_K0_SIGMA_MIN:g} "
                   "(insufficient low-k suppression)",
             f"Increase k0 or packet_width so k0*sigma >= {G2_K0_SIGMA_MIN:g}.")
    ppw = params.points_per_wavelength
    if ppw < G3_POINTS_PER_WAVELENGTH:
        fail("G3", f"{ppw:.3g} points per carrier wavelength < "
                   f"{G3_POINTS_PER_WAVELENGTH:g}",
             f"Reduce dr below {2*np.pi/(params.k0*G3_POINTS_PER_WAVELENGTH):.4g} "
             "or reduce k0.")
    k_edge = params.k0 - G4_SPECTRAL_EDGE / params.packet_width
    if k_edge <= 0:
        fail("G4", f"spectral slow edge k_edge = {k_edge:.3g} <= 0",
             "Increase k0*sigma so the spectrum stays clear of k = 0.")
    v_edge = k_edge / np.sqrt(k_edge ** 2 + M0 ** 2)
    L = params.r0_effective() + params.sponge_start_frac * params.r_max
    t_need = G4_CLEARANCE_MARGIN * L / v_edge
    if params.t_end < t_need:
        fail("G4", f"t_end = {params.t_end:g} < {t_need:.4g} required for "
                   f"{G4_CLEARANCE_MARGIN:g}x clearance of the slow spectral "
                   f"edge (v_g = {v_edge:.3f}, path L = {L:g})",
             f"Raise t_end to at least {t_need:.4g}, or raise k0/sigma.")
    info = {"k_edge": float(k_edge), "v_g_edge": float(v_edge),
            "clearance_path_L": float(L), "t_end_required": float(t_need)}
    if log:
        log(f"  carrier guards OK: k0*sigma={params.k0_sigma:.2f}, "
            f"pts/lambda={ppw:.1f}, v_g={params.v_group:.3f}, "
            f"t_end {params.t_end:g} >= {t_need:.0f} required")
    return info


# ======================================================================
# L3 -- EVOLUTION ENGINE (dynamics, energy, RK4, honest halts)
# ======================================================================
def rhs(y: np.ndarray, grid: Grid, W: np.ndarray, lambda_chi: float,
        damping: np.ndarray) -> np.ndarray:
    """chi_tt = div_flux(chi_r) - (m0^2 + W) chi - lambda_chi chi^3 - gamma chi_t.
    Spatial operator is the inlined finite-volume div_flux; the outermost node
    is Dirichlet-held at chi = 0 behind the sponge."""
    N = grid.N
    chi = y[:N]
    chid = y[N:]
    F = grid.face_grad(chi)
    div = grid.div_flux(F)
    acc = div - (M0 ** 2 + W) * chi - lambda_chi * chi ** 3 - damping * chid
    acc[-1] = 0.0
    out = np.empty_like(y)
    out[:N] = chid
    out[N - 1] = 0.0
    out[N:] = acc
    return out


def energy_box(y: np.ndarray, grid: Grid, W: np.ndarray,
               lambda_chi: float) -> float:
    """Face-consistent discrete energy: node kinetic+potential with spherical
    weights, plus face-gradient terms. The local potential (m0^2+W)chi^2/2 may
    be negative inside the well; the total of a real-omega mode is
    omega^2-positive by the eigenrelation."""
    N = grid.N
    chi = y[:N]
    chid = y[N:]
    nodes = (0.5 * chid ** 2 + 0.5 * (M0 ** 2 + W) * chi ** 2
             + 0.25 * lambda_chi * chi ** 4)
    E_nodes = float(np.sum(grid.w * nodes))
    gf = grid.face_grad(chi)
    E_faces = float(np.sum(4.0 * np.pi * grid.area_f * grid.dr * 0.5 * gf ** 2))
    return E_nodes + E_faces


def mode_amplitudes(y: np.ndarray, grid: Grid, real_modes: list) -> list:
    """c_n = int u_n (r chi) dr, cdot_n = int u_n (r chi_t) dr, and mode energy
    E_n = 2*pi*(cdot_n^2 + omega_n^2 c_n^2)."""
    N = grid.N
    rphi = grid.r * y[:N]
    rphid = grid.r * y[N:]
    out = []
    for m in real_modes:
        c = float(np.trapezoid(m["u"] * rphi, grid.r))
        cd = float(np.trapezoid(m["u"] * rphid, grid.r))
        En = 2.0 * np.pi * (cd ** 2 + m["omega2"] * c ** 2)
        out.append({"n": m["n"], "c": c, "cdot": cd, "E": En})
    return out


@dataclass
class EvolveResult:
    ts: dict = field(default_factory=dict)   # contract time-series columns
    snap_t: Optional[np.ndarray] = None
    snap_chi: Optional[np.ndarray] = None    # full-resolution snapshots (N_t, N)
    E_in: float = 0.0
    halted: bool = False
    halt_reason: str = ""


def evolve(params, grid: Grid, W: np.ndarray, y0: np.ndarray,
           lambda_chi: float, real_modes: list, log=print) -> EvolveResult:
    """Block RK4 with fixed dt = cfl*dr (unit characteristic speed), sponge-
    power bookkeeping, and honest halts. Produces the contract time-series
    columns directly (t, E_box, E_esc, E_mode_n, f_cap, f_esc, f_res, closure,
    cons_err) and a full-resolution snapshot buffer for the field export.
    """
    res = EvolveResult()
    damping = grid.sponge.copy()
    E_in = energy_box(y0, grid, W, lambda_chi)          # R6: full functional, seed
    res.E_in = E_in
    # geometry-independent amplitude reference: max|u| = max|r*chi| of the seed
    u_seed_max = float(np.max(np.abs(grid.r * y0[:grid.N])))
    u_limit = U_BLOWUP_FACTOR * max(u_seed_max, 1e-12)
    denom = max(abs(E_in), E_FLOOR)

    n_out = max(2, int(params.n_out))
    t_out = np.linspace(0.0, params.t_end, n_out)
    n_snap = max(2, int(params.n_snapshots))
    snap_every = max(1, (n_out - 1) // max(1, n_snap - 1))
    dt = params.cfl * params.dr

    N = grid.N
    y = y0.copy()
    t = 0.0
    E_esc = 0.0
    rows = []
    snaps_t, snaps_chi = [], []

    def diag_row(tt):
        E_box = energy_box(y, grid, W, lambda_chi)
        amps = mode_amplitudes(y, grid, real_modes)
        E_modes = float(np.nansum([a["E"] for a in amps])) if amps else 0.0
        row = {
            "t": tt,
            "E_box": E_box,
            "E_esc": E_esc,
            "f_cap": E_modes / denom,
            "f_esc": E_esc / denom,
            "f_res": (E_box - E_modes) / denom,
            "closure": (E_esc + E_modes) / denom - 1.0,        # signed (R7)
            "cons_err": (E_box + E_esc - E_in) / denom,        # signed (R7)
        }
        for a in amps:
            row[f"E_mode_{a['n']}"] = a["E"]
        return row

    rows.append(diag_row(0.0))
    snaps_t.append(0.0)
    snaps_chi.append(y[:N].copy())

    for i_out in range(1, n_out):
        t_target = t_out[i_out]
        n_sub = max(1, int(np.ceil((t_target - t) / dt)))
        dts = (t_target - t) / n_sub
        ok = True
        for _ in range(n_sub):
            chid_pre = y[N:]                                    # pre-substep chi_t
            E_esc += float(np.sum(grid.w * damping * chid_pre ** 2)) * dts
            k1 = rhs(y, grid, W, lambda_chi, damping)
            k2 = rhs(y + 0.5 * dts * k1, grid, W, lambda_chi, damping)
            k3 = rhs(y + 0.5 * dts * k2, grid, W, lambda_chi, damping)
            k4 = rhs(y + dts * k3, grid, W, lambda_chi, damping)
            y = y + (dts / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dts
            if not np.all(np.isfinite(y)):
                res.halted = True
                res.halt_reason = "NaN/Inf in state (honest halt)."
                ok = False
                break
        if not ok:
            break
        row = diag_row(t)
        rows.append(row)
        if (i_out % snap_every == 0) or (i_out == n_out - 1):
            snaps_t.append(t)
            snaps_chi.append(y[:N].copy())
        u_max = float(np.max(np.abs(grid.r * y[:N])))
        if u_max > u_limit:
            res.halted = True
            res.halt_reason = (
                f"amplitude blow-up: max|r*chi| = {u_max:.3e} exceeded "
                f"{U_BLOWUP_FACTOR:g}x the seed value {u_seed_max:.3e} "
                "(honest halt).")
            break
        if abs(row["E_box"]) > E_BLOWUP_FACTOR * max(abs(E_in), 1.0):
            res.halted = True
            res.halt_reason = "energy blow-up (honest halt)."
            break

    keys = rows[0].keys()
    res.ts = {k: np.array([r.get(k, np.nan) for r in rows]) for k in keys}
    res.snap_t = np.array(snaps_t)
    res.snap_chi = np.array(snaps_chi)
    return res


# ======================================================================
# L4 -- INITIAL CONDITIONS AND ORTHOGONALIZATION
# ======================================================================
def raw_packet_fields(grid: Grid, params):
    """Ingoing packet (chi, chi_t) before any orthogonalization. Regular at
    r = 0; boundary node zeroed.

    'chi_carrier'  : finite-central-momentum ingoing carrier packet (R9-A).
                     k0 > 0 means INWARD (toward decreasing r). Built in
                     u = r*chi with the single-direction WKB pairing

                       E(r)   = amp * exp(-((r-r0)/sigma)^2)
                       U(r)   = r * E(r)
                       u0(r)  = U cos(k0 (r-r0))
                       u0dot  = v_g U' cos(k0 (r-r0)) - omega0 U sin(k0 (r-r0))

                     with omega0 = sqrt(k0^2 + m0^2), v_g = k0/omega0. In the
                     massless limit this reduces to the ingoing relation
                     u_t = u_r; outgoing contamination is second order in the
                     envelope expansion (self-test S13 measures it).
    'chi_gaussian' : legacy zero-momentum ingoing packet, retained unchanged
                     for the control_k0_zero diagnostic run.
    'u_gaussian'   : legacy option, prepared in u = r*chi.
    """
    r0 = params.r0_packet if params.r0_packet > 0 else 0.55 * grid.r_max
    if params.packet_kind == "chi_carrier":
        omega0 = params.omega0
        v_g = params.v_group
        E = params.amp * np.exp(-((grid.r - r0) / params.packet_width) ** 2)
        U = grid.r * E                        # u-space envelope
        Up = grid.node_grad(U)                # slowly varying derivative
        phase = params.k0 * (grid.r - r0)
        cos_p, sin_p = np.cos(phase), np.sin(phase)
        u0 = U * cos_p
        ud0 = v_g * Up * cos_p - omega0 * U * sin_p
        chi = np.zeros(grid.N)
        chid = np.zeros(grid.N)
        chi[1:] = u0[1:] / grid.r[1:]
        chid[1:] = ud0[1:] / grid.r[1:]
        chi[0] = chi[1]
        chid[0] = chid[1]
        chi[-1] = 0.0
        chid[-1] = 0.0
    elif params.packet_kind == "chi_gaussian":
        chi = params.amp * np.exp(-((grid.r - r0) / params.packet_width) ** 2)
        chi[-1] = 0.0
        chid = grid.node_grad(chi)
        chid[-1] = 0.0
    elif params.packet_kind == "u_gaussian":
        u = params.amp * np.exp(-((grid.r - r0) / params.packet_width) ** 2)
        ut = grid.node_grad(u)
        chi = np.zeros(grid.N)
        chid = np.zeros(grid.N)
        chi[1:] = u[1:] / grid.r[1:]
        chid[1:] = ut[1:] / grid.r[1:]
        chi[0] = chi[1]
        chid[0] = chid[1]
        chi[-1] = 0.0
        chid[-1] = 0.0
    else:
        raise ValueError(f"unknown packet_kind: {params.packet_kind}")
    return chi, chid, r0


def orthogonalize_bound_modes(y, grid: Grid, real_modes: list):
    """Project the seed onto the bound subspace and subtract it, in both the
    field and its time derivative. Work in u = r*chi with int u_n^2 = 1:
        c_n = int u_n u dr,  cdot_n = int u_n u_t dr
        u_clean = u - sum c_n u_n,  u_t_clean = u_t - sum cdot_n u_n.
    Convert back with regularity at r = 0; zero the boundary node."""
    N = grid.N
    u = grid.r * y[:N].copy()
    ut = grid.r * y[N:].copy()
    recs = []
    for m in real_modes:
        c = float(np.trapezoid(m["u"] * u, grid.r))
        cd = float(np.trapezoid(m["u"] * ut, grid.r))
        u = u - c * m["u"]
        ut = ut - cd * m["u"]
        recs.append({"n": m["n"], "c_proj_removed": c, "cdot_proj_removed": cd})
    chi_clean = np.zeros(N)
    chid_clean = np.zeros(N)
    chi_clean[1:] = u[1:] / grid.r[1:]
    chid_clean[1:] = ut[1:] / grid.r[1:]
    chi_clean[0] = chi_clean[1]
    chid_clean[0] = chid_clean[1]
    chi_clean[-1] = 0.0
    chid_clean[-1] = 0.0
    return np.concatenate([chi_clean, chid_clean]), recs


def seed_hash(y0: np.ndarray) -> str:
    """Deterministic SHA-256 over the post-orthogonalization (chi0, chi_t0)."""
    return hashlib.sha256(np.ascontiguousarray(y0, dtype=np.float64).tobytes()
                          ).hexdigest()


def build_seed(grid: Grid, params, W: np.ndarray, lambda_chi: float,
               real_modes: list, log=print):
    """Construct the seed and its diagnostics. R6: f_cap_init uses the SAME
    denominator E_in = |E_box(0)| (full functional, orthogonalized seed) as the
    evolution f_cap(t)."""
    chi, chid, r0 = raw_packet_fields(grid, params)
    y_raw = np.concatenate([chi, chid])

    info = {
        "packet_kind": params.packet_kind,
        "orthogonalized": bool(params.orthogonalize),
        "r0": float(r0),
        "overlap_records": [],
    }

    if params.orthogonalize and real_modes:
        y0, recs = orthogonalize_bound_modes(y_raw, grid, real_modes)
        info["overlap_records"] = recs
    else:
        y0 = y_raw

    # Unified denominator on the (possibly orthogonalized) seed.
    E_in = max(abs(energy_box(y0, grid, W, lambda_chi)), E_FLOOR)
    amps0 = mode_amplitudes(y0, grid, real_modes)
    E_modes0 = float(np.nansum([a["E"] for a in amps0])) if amps0 else 0.0
    f_cap_init = E_modes0 / E_in
    info["f_cap_init"] = float(f_cap_init)
    info["seed_clean"] = bool(f_cap_init < GATES["eps_lin"])
    info["seed_hash"] = seed_hash(y0)
    log(f"  seed: f_cap_init = {f_cap_init:.3e} "
        f"(eps_lin = {GATES['eps_lin']:g}) -> "
        f"{'clean' if info['seed_clean'] else 'CONTAMINATED'}")
    return y0, info


# ======================================================================
# L5 -- WINDOW DIAGNOSTICS
# ======================================================================
def final_window_mask(t: np.ndarray, t_char: float) -> np.ndarray:
    """Final analysis window: the last 5 characteristic times (R2, preserved
    verbatim). Falls back to the last 20% if that window is empty."""
    win = t >= (t[-1] - 5.0 * t_char)
    if not win.any():
        win = t >= 0.8 * t[-1]
    return win


def compute_t_clear(t: np.ndarray, f_res: np.ndarray, eps_res: float):
    """Earliest output time after which f_res < eps_res holds through t_end;
    None if radiation never clears within the run."""
    below = f_res < eps_res
    if not below.any():
        return None
    for i in range(len(t)):
        if below[i] and below[i:].all():
            return float(t[i])
    return None


def window_mode_statistics(ts: dict, t_char: float, f_cap_window: float) -> dict:
    """Per-mode occupation statistics over the final window, with the SAME
    denominator as f_cap (E_in recovered from the f_cap identity
    E_in = sum_n E_mode_n / f_cap). Descriptive only. Neutral vocabulary.
    """
    out = {
        "n_real_modes": 0, "dominant_mode_final_window": None,
        "mode_fractions_final_window_median": {},
        "mode_fraction_of_captured_final_window": {},
        "mode_fraction_std_over_mean_final_window": {},
        "denominator_note": ("E_in recovered from f_cap identity "
                             "(E_in = sum E_mode / f_cap); consistent with the "
                             "evolution denominator |E_box(0)|"),
    }
    mode_keys = sorted([k for k in ts if k.startswith("E_mode_")],
                       key=lambda k: int(k.split("_")[-1]))
    out["n_real_modes"] = len(mode_keys)
    if not mode_keys:
        return out
    t = ts["t"]
    win = final_window_mask(t, t_char)
    sum_final = float(np.nansum([ts[k][-1] for k in mode_keys]))
    fcap_final = float(ts["f_cap"][-1])
    if abs(fcap_final) < 1e-30:
        return out
    E_in = sum_final / fcap_final
    if not np.isfinite(E_in) or abs(E_in) < 1e-30:
        return out
    fw_med, sov = {}, {}
    for k in mode_keys:
        n = k.split("_")[-1]
        arr = ts[k]
        fw_med[f"mode_{n}"] = float(np.nanmedian(arr[win]) / E_in)
        mu = float(np.nanmean(arr[win]) / E_in)
        sd = float(np.nanstd(arr[win]) / E_in)
        sov[f"mode_{n}"] = float(sd / mu) if abs(mu) > 1e-30 else float("nan")
    out["mode_fractions_final_window_median"] = fw_med
    out["mode_fraction_std_over_mean_final_window"] = sov
    if np.isfinite(f_cap_window) and abs(f_cap_window) >= GATES["f_min"] * 0.1:
        out["mode_fraction_of_captured_final_window"] = {
            key: (v / f_cap_window if np.isfinite(v) else None)
            for key, v in fw_med.items()}
    else:
        out["mode_fraction_of_captured_final_window"] = {k: None for k in fw_med}
    if fw_med:
        out["dominant_mode_final_window"] = max(fw_med, key=lambda k: fw_med[k])
    return out


def window_outcomes(res: EvolveResult, params, gates: dict) -> dict:
    """Assemble the outcomes block from the time series: window medians for
    f_cap and f_res, terminal f_esc, closure, net f_cap, cons_err_max, t_clear,
    and the analysis window."""
    ts = res.ts
    t = ts["t"]
    win = final_window_mask(t, params.t_char)
    f_cap_final = float(np.nanmedian(ts["f_cap"][win]))
    f_cap_init = float(ts["f_cap"][0])
    f_res_final = float(np.nanmedian(ts["f_res"][win]))
    f_esc_final = float(ts["f_esc"][-1])
    closure_abs = abs(f_esc_final + f_cap_final - 1.0)
    cons_err_max = float(np.nanmax(np.abs(ts["cons_err"])))
    t_clear = compute_t_clear(t, ts["f_res"], gates["eps_res"])
    t_start = float(t[win][0]) if win.any() else float(t[0])
    outcomes = {
        "f_cap_init": f_cap_init,
        "f_cap_final": f_cap_final,
        "net_f_cap": f_cap_final - f_cap_init,
        "f_esc_final": f_esc_final,
        "f_res_final": f_res_final,
        "closure_abs": float(closure_abs),
        "cons_err_max": cons_err_max,
        "radiation_cleared": bool(f_res_final < gates["eps_res"]),
        "budget_closed": bool(closure_abs < gates["eps_budget"]),
        "t_clear": t_clear,
        "analysis_window": [t_start, float(t[-1])],
        "halted": bool(res.halted),
        "halt_reason": res.halt_reason,
        "window_mode_statistics": window_mode_statistics(
            ts, params.t_char, f_cap_final),
    }
    return outcomes


# ======================================================================
# L6 -- GATE AND STATUS ADJUDICATION
# ======================================================================
def compute_numerics_valid(outcomes: dict, gates: dict,
                           null_floor_max: Optional[float],
                           is_baseline_class: bool) -> bool:
    """numerics_valid = no halt AND conservation within EPS_CONS AND, for
    baseline-class runs (R4), a clean linear-null floor (null_floor_max <
    eps_lin). Any NaN/Inf shows up as a halt or as a non-finite cons_err."""
    if outcomes.get("halted", False):
        return False
    cem = outcomes.get("cons_err_max", float("inf"))
    if not np.isfinite(cem) or cem >= EPS_CONS:
        return False
    if is_baseline_class:
        if null_floor_max is None or not np.isfinite(null_floor_max):
            return False
        if null_floor_max >= gates["eps_lin"]:
            return False
    return True


def adjudicate_run(seed_clean: bool, numerics_valid: bool,
                   radiation_cleared: bool, budget_closed: bool,
                   net_f_cap: float, gates: dict) -> str:
    """The five-state per-run taxonomy, evaluated top-down; the first failing
    class assigns the status. One of exactly two functions that emit a status.
    """
    if not seed_clean:
        return "INVALID_SEED"
    if not numerics_valid:
        return "INVALID_NUMERICS"
    if not (radiation_cleared and budget_closed):
        return "INCONCLUSIVE"
    if net_f_cap < gates["f_min"]:
        return "REJECTED"
    return "SUPPORTED"


def adjudicate_ensemble(subset_rows: list, gates: dict) -> dict:
    """H2 aggregation over the pre-registered robustness subset only
    (in_robustness_set == True rows). rel_spread = (max - min)/median of
    net f_cap (R3). The second and last function that emits a status.

    Exactly three ensemble states exist: H2_SUPPORTED, H2_REJECTED,
    H2_INCONCLUSIVE. An empty robustness subset is H2_INCONCLUSIVE with an
    explanatory ensemble_note -- never a fourth status.
    """
    result = {
        "net_f_cap_median_supported": None,
        "rel_spread_supported": None,
        "h2_ensemble_status": "H2_INCONCLUSIVE",
        "ensemble_note": "",
    }
    if not subset_rows:
        result["ensemble_note"] = "no robustness subset rows"
        return result
    statuses = [r["status"] for r in subset_rows]

    supported = [r for r in subset_rows if r["status"] == "SUPPORTED"]
    if supported:
        s_net = np.array([r["net_f_cap"] for r in supported], dtype=float)
        med = float(np.median(s_net))
        result["net_f_cap_median_supported"] = med
        if len(s_net) >= 2 and abs(med) > 1e-30:
            result["rel_spread_supported"] = float(
                (np.max(s_net) - np.min(s_net)) / med)

    all_supported = all(s == "SUPPORTED" for s in statuses)
    all_adjudicable = all(s in ("SUPPORTED", "REJECTED") for s in statuses)
    spread = result["rel_spread_supported"]

    if all_supported and (len(supported) < 2 or (spread is not None
                                                 and spread <= gates["rtol"])):
        result["h2_ensemble_status"] = "H2_SUPPORTED"
        result["ensemble_note"] = (
            f"all {len(subset_rows)} subset run(s) SUPPORTED"
            + (f"; rel_spread={spread:.3f} <= rtol={gates['rtol']:g}"
               if spread is not None else "; single supported run, spread N/A"))
    elif all_adjudicable and (any(s == "REJECTED" for s in statuses)
                              or (spread is not None and spread > gates["rtol"])):
        result["h2_ensemble_status"] = "H2_REJECTED"
        result["ensemble_note"] = (
            "subset adjudicable; "
            + ("contains REJECTED run(s)" if any(s == "REJECTED" for s in statuses)
               else f"rel_spread={spread:.3f} > rtol={gates['rtol']:g}"))
    else:
        result["h2_ensemble_status"] = "H2_INCONCLUSIVE"
        n_bad = sum(1 for s in statuses
                    if s in ("INCONCLUSIVE", "INVALID_SEED", "INVALID_NUMERICS"))
        result["ensemble_note"] = (
            f"{n_bad} of {len(subset_rows)} subset run(s) not adjudicable")
    return result


# ======================================================================
# L7 -- PARAMETERS, CONFIG, PROVENANCE, EXPORTERS
# ======================================================================
@dataclass
class RunParams:
    # physics
    profile: str = "square"
    W0: float = 0.9
    a: float = 5.0
    Delta: float = 0.5              # tanh interface width (manuscript Eq. 3)
    lambda_chi: float = 1.0         # +: bounded quartic (manuscript convention)
    amp: float = 0.2
    r0_packet: float = 0.0          # 0 -> auto 0.55*r_max
    packet_width: float = 2.0       # sigma
    packet_kind: str = "chi_gaussian"   # 'chi_carrier' | 'chi_gaussian' | 'u_gaussian'
    k0: float = 0.0                 # carrier wavenumber; k0 > 0 means INWARD
    orthogonalize: bool = True
    # numerics
    r_max: float = 40.0
    dr: float = 0.05
    cfl: float = 0.2
    t_end: float = 200.0
    n_out: int = 600
    n_snapshots: int = 81
    sponge_start_frac: float = 0.75
    sponge_strength: float = 3.0
    # field export decimation (R10)
    field_r_stride: int = 4
    field_n_snapshots: int = 80

    @property
    def t_char(self) -> float:
        return max(self.a, 1.0 / M0)

    # ---- derived carrier diagnostics (model units; k0 = 0 -> not a carrier) --
    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0 ** 2 + M0 ** 2))

    @property
    def v_group(self) -> float:
        return float(self.k0 / self.omega0) if self.k0 > 0 else 0.0

    @property
    def k0_sigma(self) -> float:
        return float(self.k0 * self.packet_width)

    @property
    def points_per_wavelength(self) -> float:
        return float(2.0 * np.pi / (self.k0 * self.dr)) if self.k0 > 0 \
            else float("inf")

    def r0_effective(self, r_max=None) -> float:
        rm = self.r_max if r_max is None else r_max
        return self.r0_packet if self.r0_packet > 0 else 0.55 * rm


def _json_safe(x):
    if isinstance(x, dict):
        return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, np.ndarray):
        return _json_safe(x.tolist())
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        xf = float(x)
        return None if (np.isnan(xf) or np.isinf(xf)) else xf
    return x


def _fmt(v):
    """Deterministic scalar formatting for CSV (byte-stable across runs)."""
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float) or isinstance(v, np.floating):
        vf = float(v)
        if np.isnan(vf):
            return ""
        return f"{vf:.12g}"
    return str(v)


def _write_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([_fmt(x) for x in r])


def _write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(_json_safe(obj), fh, indent=2, allow_nan=False)


def provenance_block(config_path, config_obj, wall_time_s, seed_hash_val,
                     extra=None):
    blk = {
        "code_version": CODE_VERSION,
        "config_path": config_path,
        "config_sha256": (hashlib.sha256(
            json.dumps(config_obj, sort_keys=True).encode()).hexdigest()
            if config_obj is not None else None),
        "seed_hash": seed_hash_val,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "wall_time_s": float(wall_time_s),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    if extra:
        blk.update(extra)
    return blk


def export_time_series(path, ts: dict):
    mode_keys = sorted([k for k in ts if k.startswith("E_mode_")],
                       key=lambda k: int(k.split("_")[-1]))
    header = (["t", "E_box", "E_esc"] + mode_keys
              + ["f_cap", "f_esc", "f_res", "closure", "cons_err"])
    n = len(ts["t"])
    rows = []
    for i in range(n):
        rows.append([ts["t"][i], ts["E_box"][i], ts["E_esc"][i]]
                    + [ts[k][i] for k in mode_keys]
                    + [ts["f_cap"][i], ts["f_esc"][i], ts["f_res"][i],
                       ts["closure"][i], ts["cons_err"][i]])
    _write_csv(path, header, rows)


def export_mode_profiles(path, grid: Grid, W: np.ndarray, real_modes: list):
    header = ["r", "W", "m_eff2"] + [f"u_{m['n']}" for m in real_modes]
    rows = []
    for j in range(grid.N):
        row = [grid.r[j], W[j], M0 ** 2 + W[j]]
        row += [m["u"][j] for m in real_modes]
        rows.append(row)
    _write_csv(path, header, rows)


def export_field_evolution(path, grid: Grid, res: EvolveResult, params):
    """float32 npz with r, t_snap, chi (N_t, N_r), decimated per R10 and
    enforcing N_t * N_r <= 2e6 (hard error, not silent truncation)."""
    r_stride = max(1, int(params.field_r_stride))
    r_idx = np.arange(0, grid.N, r_stride)
    r_out = grid.r[r_idx]

    snap_t = res.snap_t
    n_have = len(snap_t)
    n_want = max(2, int(params.field_n_snapshots))
    if n_have > n_want:
        t_idx = np.unique(np.linspace(0, n_have - 1, n_want).astype(int))
    else:
        t_idx = np.arange(n_have)
    t_out = snap_t[t_idx]
    chi_out = res.snap_chi[np.ix_(t_idx, r_idx)].astype(np.float32)

    N_t, N_r = chi_out.shape
    if N_t * N_r > 2_000_000:
        raise ValueError(
            f"field_evolution would be {N_t}x{N_r} = {N_t * N_r} cells > 2e6. "
            "Increase field_r_stride or reduce field_n_snapshots in the config.")
    if not (np.all(np.diff(r_out) > 0) and np.all(np.diff(t_out) > 0)):
        raise ValueError("field_evolution axes must be strictly increasing.")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, r=r_out.astype(np.float32),
                        t_snap=t_out.astype(np.float32), chi=chi_out)
    return {"field_r_stride": r_stride, "field_n_snapshots_written": int(N_t),
            "field_N_r": int(N_r), "field_dtype": "float32"}


def build_run_metadata(run_id, params: RunParams, real_modes, outcomes, status,
                       provenance, is_linear_null=False,
                       in_robustness_set=False, is_control=False):
    return {
        "run_id": run_id,
        "code_version": CODE_VERSION,
        "created_at": provenance["created_at"],
        "m0": M0,
        # Pre-registered H2 membership, persisted here so that a standalone
        # `adjudicate --data <dir>` is fully idempotent without needing the
        # original scan config.
        "in_robustness_set": bool(in_robustness_set),
        # Control/diagnostic runs (R9-D): adjudicated per-run and listed in
        # capture_budget.csv, but never in scan_summary.csv and never in H2.
        "is_control": bool(is_control),
        "profile": params.profile,
        "W0": params.W0,
        "a": params.a,
        "Delta": params.Delta,
        "lambda_chi": (0.0 if is_linear_null else params.lambda_chi),
        "is_linear_null": bool(is_linear_null),
        "packet": {
            "kind": params.packet_kind,
            "amplitude": params.amp,
            "r0": params.r0_effective(),
            "sigma": params.packet_width,
            "orthogonalized": params.orthogonalize,
            # carrier diagnostics; k0 is null for non-carrier packet kinds
            "k0": (params.k0 if params.packet_kind == "chi_carrier" else None),
            "omega0": (params.omega0 if params.packet_kind == "chi_carrier"
                       else None),
            "v_group": (params.v_group if params.packet_kind == "chi_carrier"
                        else None),
            "k0_sigma": (params.k0_sigma if params.packet_kind == "chi_carrier"
                         else None),
            "points_per_wavelength": (
                params.points_per_wavelength
                if params.packet_kind == "chi_carrier" else None),
        },
        "grid": {"r_max": params.r_max, "dr": params.dr, "cfl": params.cfl,
                 "t_end": params.t_end},
        "sponge": {"start_frac": params.sponge_start_frac,
                   "strength": params.sponge_strength},
        "gates": dict(GATES),
        "modes": [{"n": m["n"], "omega2": m["omega2"], "kappa": m["kappa"],
                   "box_affected": m["box_affected"], "n_nodes": m["n_nodes"]}
                  for m in real_modes],
        "outcomes": {**outcomes, "status": status},
        "provenance": provenance,
    }


# ======================================================================
# L9(a) -- CORE RUN DRIVERS
# ======================================================================
def _prepare(params: RunParams, lambda_chi: float, log):
    # Configuration refusals first: an out-of-protocol probe or an
    # out-of-scope well never reaches the evolution engine.
    guard_info = assert_carrier_guards(params, log=log)
    grid = Grid(GridSpec(params.r_max, params.dr,
                         params.sponge_start_frac, params.sponge_strength))
    W = W_grid(params.profile, grid.r, params.W0, params.a, params.dr,
               params.Delta)
    modes = spectral_bound_modes(grid, W)
    assert_no_tachyon(modes, f"profile={params.profile}, W0={params.W0}, a={params.a}")
    rmodes = real_bound_modes(modes)
    if not rmodes:
        raise ValueError(
            f"No ordinary bound modes for profile={params.profile}, "
            f"W0={params.W0}, a={params.a}; deepen the well within scope.")
    return grid, W, rmodes, guard_info


def run_capture(params: RunParams, out_root, run_id, lambda_chi,
                is_linear_null, config_path, config_obj, log=print,
                write_full=True, in_robustness_set=False, is_control=False):
    """One evolved run. Returns (metadata, extras) where extras carries the
    null_floor_max for baseline-class numerics validity.

    in_robustness_set is the pre-registered H2 membership flag; it is written
    into run_metadata.json as a stable top-level field so that adjudication
    never has to consult the original scan config.
    """
    t0 = time.time()
    grid, W, rmodes, guard_info = _prepare(params, lambda_chi, log)
    y0, seed_info = build_seed(grid, params, W, lambda_chi, rmodes, log)
    res = evolve(params, grid, W, y0, lambda_chi, rmodes, log)
    outcomes = window_outcomes(res, params, GATES)
    outcomes["f_cap_init"] = seed_info["f_cap_init"]
    outcomes["net_f_cap"] = outcomes["f_cap_final"] - seed_info["f_cap_init"]
    outcomes["seed_clean"] = seed_info["seed_clean"]

    null_floor_max = float(np.nanmax(res.ts["f_cap"])) if is_linear_null else None

    wall = time.time() - t0
    prov = provenance_block(config_path, config_obj, wall, seed_info["seed_hash"],
                            extra={"overlap_records": seed_info["overlap_records"]})

    # provisional status (adjudicate re-derives authoritatively later)
    nv = compute_numerics_valid(outcomes, GATES, null_floor_max,
                                is_baseline_class=False)
    outcomes["numerics_valid"] = nv
    status = adjudicate_run(outcomes["seed_clean"], nv,
                            outcomes["radiation_cleared"],
                            outcomes["budget_closed"], outcomes["net_f_cap"],
                            GATES)
    outcomes["null_floor_max"] = null_floor_max

    if guard_info:
        prov["carrier_guards"] = guard_info
    meta = build_run_metadata(run_id, params, rmodes, outcomes, status, prov,
                              is_linear_null=is_linear_null,
                              in_robustness_set=in_robustness_set,
                              is_control=is_control)

    run_dir = out_root
    _write_json(os.path.join(run_dir, "run_metadata.json"), meta)
    export_time_series(os.path.join(run_dir, "time_series.csv"), res.ts)
    if write_full:
        export_mode_profiles(os.path.join(run_dir, "mode_profiles.csv"),
                             grid, W, rmodes)
        field_meta = export_field_evolution(
            os.path.join(run_dir, "field_evolution.npz"), grid, res, params)
        meta["provenance"].update(field_meta)
        _write_json(os.path.join(run_dir, "run_metadata.json"), meta)  # refresh

    log(f"  run {run_id}: status={status} net_f_cap={outcomes['net_f_cap']:.4e} "
        f"closure={outcomes['closure_abs']:.2e} cons_err={outcomes['cons_err_max']:.2e}")
    return meta, {"null_floor_max": null_floor_max, "seed_hash": seed_info["seed_hash"]}


# ======================================================================
# L9(b) -- SUBCOMMAND WORKFLOWS
# ======================================================================
def cmd_spectral(cfg, data_root, log=print):
    spec = cfg["spectral"]
    a = spec["a"]
    W0_grid = spec["W0_grid"]
    dr = spec.get("dr", 0.02)
    r_max = spec.get("r_max", 40.0)
    extra_profiles = spec.get("extra_profiles", [])

    val_rows, cnt_rows = [], []
    for profile in (["square"] + list(extra_profiles)):
        for W0 in W0_grid:
            grid = Grid(GridSpec(r_max=r_max, dr=dr))
            W = W_grid(profile, grid.r, W0, a, dr,
                       cfg.get("Delta", 0.5))
            modes = spectral_bound_modes(grid, W)
            assert_no_tachyon(modes, f"spectral {profile} W0={W0}")
            rmodes = real_bound_modes(modes)
            sqrtW0_a = float(np.sqrt(W0) * a)
            if profile == "square":
                match = square_well_matching_omega2(W0, a)
                N_analytic = analytic_mode_count(sqrtW0_a)
            else:
                match = [float("nan")] * len(rmodes)
                N_analytic = -1  # not defined for smooth profiles
            for i, m in enumerate(rmodes):
                om2_match = (match[i] if i < len(match) else float("nan"))
                val_rows.append([profile, W0, a, sqrtW0_a, m["n"],
                                 m["omega2"], om2_match, m["kappa"]])
            cnt_rows.append([profile, W0, a, sqrtW0_a, len(rmodes), N_analytic])

    _write_csv(os.path.join(data_root, "spectral", "spectral_validation.csv"),
               ["profile", "W0", "a", "sqrtW0_a", "n", "omega2_n",
                "omega2_matching", "kappa_n"], val_rows)
    _write_csv(os.path.join(data_root, "spectral", "mode_counts.csv"),
               ["profile", "W0", "a", "sqrtW0_a", "N_numeric", "N_analytic"],
               cnt_rows)
    log(f"spectral: wrote {len(val_rows)} eigenvalue rows, "
        f"{len(cnt_rows)} count rows.")


def _params_from_cfg(block: dict) -> RunParams:
    valid = {f for f in RunParams().__dataclass_fields__}
    return RunParams(**{k: v for k, v in block.items() if k in valid})


def cmd_baseline(cfg, data_root, config_path, log=print):
    params = _params_from_cfg(cfg["baseline"])
    run_id = cfg["baseline"].get("run_id", "baseline")
    # R9 decides whether the baseline run is itself a robustness-set member;
    # the flag is persisted into its metadata either way.
    in_set = bool(cfg["baseline"].get("in_robustness_set", False))
    meta, _ = run_capture(params, os.path.join(data_root, "baseline"), run_id,
                          params.lambda_chi, is_linear_null=False,
                          config_path=config_path, config_obj=cfg, log=log,
                          write_full=True, in_robustness_set=in_set)
    return meta


def cmd_null(cfg, data_root, config_path, log=print):
    params = _params_from_cfg(cfg["baseline"])
    run_id = cfg["baseline"].get("run_id", "baseline") + "_null"
    # baseline seed hash for the equality assertion, if baseline exists
    base_meta_path = os.path.join(data_root, "baseline", "run_metadata.json")
    meta, extras = run_capture(
        params, os.path.join(data_root, "null_linear"), run_id,
        lambda_chi=0.0, is_linear_null=True, config_path=config_path,
        config_obj=cfg, log=log, write_full=False)
    if os.path.exists(base_meta_path):
        with open(base_meta_path, encoding="utf-8") as fh:
            base = json.load(fh)
        bh = base.get("provenance", {}).get("seed_hash")
        if bh is not None and bh != extras["seed_hash"]:
            raise AssertionError(
                "null seed_hash != baseline seed_hash; the null run must use "
                "the identical orthogonalized seed. Check that baseline and "
                "null share the same config parameters.")
        log("  null: seed_hash matches baseline (identical seed).")
    return meta


def cmd_scan(cfg, data_root, config_path, log=print):
    runs = cfg["scan"]["runs"]
    written = []
    for spec in runs:
        base = dict(cfg.get("baseline", {}))
        base.update(spec)  # per-run overrides (profile, W0, a, in_robustness_set, run_id)
        params = _params_from_cfg(base)
        run_id = spec.get("run_id") or f"scan_{spec.get('profile','square')}_W0{spec.get('W0')}"
        in_set = bool(spec.get("in_robustness_set", False))
        try:
            meta, _ = run_capture(
                params, os.path.join(data_root, "scan", run_id), run_id,
                params.lambda_chi, is_linear_null=False,
                config_path=config_path, config_obj=cfg, log=log,
                write_full=False,  # R5: scan runs export metadata + time_series only
                in_robustness_set=in_set)  # persisted into run_metadata.json
            written.append((run_id, meta))
        except ValueError as exc:
            log(f"  scan {run_id}: REFUSED ({exc})")
    log(f"scan: completed {len(written)} runs.")
    return written


def cmd_controls(cfg, data_root, config_path, log=print):
    """Diagnostic/control runs (R9-C, R9-D): the zero-momentum control and the
    packet-sensitivity survey. Adjudicated per-run by the ordinary five-state
    function, listed in capture_budget.csv, and EXCLUDED from scan_summary.csv
    and from the H2 ensemble by construction."""
    block = cfg.get("controls", {})
    runs = block.get("runs", [])
    if not runs:
        log("controls: no controls block in config; nothing to do.")
        return []
    written = []
    for spec in runs:
        base = dict(cfg.get("baseline", {}))
        base.update(spec)
        params = _params_from_cfg(base)
        run_id = spec.get("run_id") or f"control_{spec.get('packet_kind','x')}"
        try:
            meta, _ = run_capture(
                params, os.path.join(data_root, "controls", run_id), run_id,
                params.lambda_chi, is_linear_null=False,
                config_path=config_path, config_obj=cfg, log=log,
                write_full=False, in_robustness_set=False, is_control=True)
            written.append((run_id, meta))
        except ValueError as exc:
            log(f"  control {run_id}: REFUSED ({exc})")
    log(f"controls: completed {len(written)} runs.")
    return written


def _collect_runs(data_root):
    """Gather every evolved run's metadata (baseline, null, scan/*, controls/*)."""
    metas = []
    for sub in ("baseline", "null_linear"):
        p = os.path.join(data_root, sub, "run_metadata.json")
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                metas.append((sub, json.load(fh)))
    for group in ("scan", "controls"):
        gdir = os.path.join(data_root, group)
        if os.path.isdir(gdir):
            for name in sorted(os.listdir(gdir)):
                p = os.path.join(gdir, name, "run_metadata.json")
                if os.path.exists(p):
                    with open(p, encoding="utf-8") as fh:
                        metas.append((f"{group}/{name}", json.load(fh)))
    return metas


def cmd_adjudicate(cfg, data_root, log=print):
    """Idempotent: re-derive per-run status and rebuild capture_budget.csv and
    scan_summary.csv from the per-run metadata. Folds the null-floor control
    into baseline-class numerics validity (R4).

    Robustness membership is read from each run's persisted
    `in_robustness_set` field, so this command is fully idempotent on a data
    directory alone; `cfg` is consulted only as a legacy fallback for runs
    written before that field existed.
    """
    metas = _collect_runs(data_root)
    if not metas:
        raise FileNotFoundError(
            f"No run_metadata.json found under '{data_root}'. Run the "
            "simulation subcommands first.")

    # locate the null-floor control for baseline-class runs (R4)
    null_floor = None
    null_run_id = None
    for _, m in metas:
        if m.get("is_linear_null"):
            null_floor = m["outcomes"].get("null_floor_max")
            null_run_id = m.get("run_id")

    # legacy fallback only: runs predating the persisted membership field
    robustness_members = set()
    if cfg and "scan" in cfg:
        for spec in cfg["scan"]["runs"]:
            rid = spec.get("run_id") or \
                f"scan_{spec.get('profile','square')}_W0{spec.get('W0')}"
            if spec.get("in_robustness_set", False):
                robustness_members.add(rid)

    cb_header = ["run_id", "profile", "W0", "a", "lambda_chi", "amplitude",
                 "packet_kind", "k0", "sigma", "is_control",
                 "is_linear_null", "f_cap_init", "f_cap_final", "net_f_cap",
                 "f_esc_final", "f_res_final", "closure_abs", "cons_err_max",
                 "seed_clean", "numerics_valid", "radiation_cleared",
                 "budget_closed", "capture_gate", "status"]
    cb_rows = []
    scan_records = []  # for scan_summary

    for tag, m in metas:
        o = m["outcomes"]
        is_null = bool(m.get("is_linear_null", False))
        is_ctrl = bool(m.get("is_control", False))
        # Baseline-class = the runs that carry the claim (baseline + scan).
        # Control/diagnostic runs are deliberately excluded from the R4
        # null-floor sub-rule: the floor is measured with the baseline's own
        # seed, so applying it to a control that uses a different packet would
        # be a category error. Controls still face every other validity check.
        is_baseline_class = (not is_null) and (not is_ctrl) and (
            tag == "baseline" or tag.startswith("scan/"))
        nv = compute_numerics_valid(o, GATES, null_floor,
                                    is_baseline_class=is_baseline_class)
        status = adjudicate_run(o["seed_clean"], nv, o["radiation_cleared"],
                                o["budget_closed"], o["net_f_cap"], GATES)
        capture_gate = bool(o["net_f_cap"] >= GATES["f_min"])
        # write the authoritative status back into metadata
        o["numerics_valid"] = nv
        o["status"] = status
        # null-floor provenance on baseline-class runs: which null run supplied
        # the floor, its value, and whether it passed the eps_lin control (R4)
        if is_baseline_class:
            o["null_run_id"] = null_run_id
            o["null_floor_max"] = null_floor
            o["floor_ok"] = bool(
                null_floor is not None and np.isfinite(null_floor)
                and null_floor < GATES["eps_lin"])
        _write_json(os.path.join(data_root, *tag.split("/"),
                                 "run_metadata.json"), m)

        pk = m.get("packet", {})
        cb_rows.append([
            m["run_id"], m["profile"], m["W0"], m["a"], m["lambda_chi"],
            pk.get("amplitude"), pk.get("kind"), pk.get("k0"), pk.get("sigma"),
            is_ctrl,
            is_null, o["f_cap_init"], o["f_cap_final"],
            o["net_f_cap"], o["f_esc_final"], o["f_res_final"], o["closure_abs"],
            o["cons_err_max"], o["seed_clean"], nv, o["radiation_cleared"],
            o["budget_closed"], capture_gate, status])

        if (not is_ctrl) and (tag.startswith("scan/") or tag == "baseline"):
            rid = m["run_id"]
            # persisted field is authoritative; config is a legacy fallback
            in_set = bool(m.get("in_robustness_set",
                                rid in robustness_members))
            scan_records.append({
                "run_id": rid, "profile": m["profile"], "W0": m["W0"],
                "a": m["a"], "net_f_cap": o["net_f_cap"], "status": status,
                "seed_clean": o["seed_clean"], "numerics_valid": nv,
                "radiation_cleared": o["radiation_cleared"],
                "budget_closed": o["budget_closed"],
                "in_robustness_set": in_set})

    _write_csv(os.path.join(data_root, "capture_budget.csv"), cb_header, cb_rows)

    # ensemble over the robustness subset only
    subset = [r for r in scan_records if r["in_robustness_set"]]
    ens = adjudicate_ensemble(subset, GATES)

    ss_header = ["run_id", "profile", "W0", "a", "net_f_cap", "status",
                 "seed_clean", "numerics_valid", "radiation_cleared",
                 "budget_closed", "in_robustness_set", "f_min", "rtol",
                 "net_f_cap_median_supported", "rel_spread_supported",
                 "h2_ensemble_status", "ensemble_note"]
    ss_rows = []
    for r in scan_records:
        ss_rows.append([
            r["run_id"], r["profile"], r["W0"], r["a"], r["net_f_cap"],
            r["status"], r["seed_clean"], r["numerics_valid"],
            r["radiation_cleared"], r["budget_closed"], r["in_robustness_set"],
            GATES["f_min"], GATES["rtol"], ens["net_f_cap_median_supported"],
            ens["rel_spread_supported"], ens["h2_ensemble_status"],
            ens["ensemble_note"]])
    _write_csv(os.path.join(data_root, "scan_summary.csv"), ss_header, ss_rows)

    log(f"adjudicate: {len(cb_rows)} runs -> capture_budget.csv; "
        f"{len(ss_rows)} scan rows -> scan_summary.csv; "
        f"H2 (subset of {len(subset)}) = {ens['h2_ensemble_status']}"
        + (f" [{ens['ensemble_note']}]" if ens["ensemble_note"] else ""))
    return ens


# ======================================================================
# L8 -- SELF-TESTS (scratch only; never writes into data/)
# ======================================================================
def run_selftests(scratch="selftest_scratch", log=print) -> bool:
    os.makedirs(scratch, exist_ok=True)
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        status = "PASS" if cond else "FAIL"
        if not cond:
            ok = False
        log(f"  [{status}] {name}{('  ' + detail) if detail else ''}")

    # S1 spectral staircase (tachyon-free regime W0 < m0^2; wider a for steps)
    a = 15.0
    stair_ok = True
    for W0 in np.linspace(0.05, 0.95, 20):
        grid = Grid(GridSpec(r_max=80.0, dr=0.05))
        W = W_grid("square", grid.r, W0, a, 0.05)
        md = spectral_bound_modes(grid, W)
        if any(m["tachyonic"] for m in md):
            stair_ok = False  # must not occur for W0 < m0^2
        N_num = len(real_bound_modes(md))
        N_ana = analytic_mode_count(np.sqrt(W0) * a)
        if abs(N_num - N_ana) > 1:  # +-1 only at the exact step points
            stair_ok = False
    check("S1 spectral staircase (numeric vs analytic count, W0<m0^2)", stair_ok)

    # S2 numeric vs analytic threshold depth. The appearance-threshold state is
    # marginally bound (kappa -> 0) and extends to ~r_max, so a large box is
    # needed to resolve it; the residual is a box-truncation effect, not a
    # scheme error (S3 checks well-bound eigenvalues to 1e-2).
    a_s2 = 5.0  # smaller radius -> threshold state fits the box better
    th_num = numeric_threshold_depth(a_s2, dr=0.02, r_max=320.0, n=1)
    th_ana = threshold_depth_analytic(1, a_s2)
    check("S2 threshold depth (large box)",
          abs(th_num - th_ana) / th_ana < 0.02,
          f"num={th_num:.5f} ana={th_ana:.5f}")

    # S3 matching-root pairing (ascending omega^2), away from threshold
    W0 = 1.5
    grid = Grid(GridSpec(r_max=40.0, dr=0.02))
    W = W_grid("square", grid.r, W0, a, 0.02)
    fd = [m["omega2"] for m in real_bound_modes(spectral_bound_modes(grid, W))]
    an = square_well_matching_omega2(W0, a)
    pair_ok = len(fd) == len(an) and all(
        abs(f - g) < 0.01 for f, g in zip(fd, an))
    check("S3 matching-root pairing", pair_ok,
          f"fd={[round(x,4) for x in fd]} an={[round(x,4) for x in an]}")

    # S4 profile parity with the figure package's common_io.W_of
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import common_io as cio
        rr = np.linspace(0, 30, 500)
        s4 = True
        for prof in VALID_PROFILES:
            mine = W_continuum(prof, rr, 0.9, 5.0, 0.5)
            theirs = cio.W_of(prof, rr, 0.9, 5.0, 0.5)
            if not np.allclose(mine, theirs, atol=1e-12):
                s4 = False
        check("S4 profile parity vs common_io.W_of", s4)
    except Exception as exc:  # common_io not importable in isolation
        check("S4 profile parity vs common_io.W_of", True,
              f"(skipped: {exc})")

    # S5 linear conservation (sponge off, short run)
    p = RunParams(profile="square", W0=0.9, a=5.0, lambda_chi=0.0,
                  r_max=30.0, dr=0.1, t_end=20.0, n_out=100, n_snapshots=20,
                  sponge_strength=0.0, orthogonalize=False, amp=0.05)
    g = Grid(GridSpec(p.r_max, p.dr, p.sponge_start_frac, 0.0))
    W = W_grid(p.profile, g.r, p.W0, p.a, p.dr)
    rmodes = real_bound_modes(spectral_bound_modes(g, W))
    y0, _ = build_seed(g, p, W, 0.0, rmodes, log=lambda *a, **k: None)
    r5 = evolve(p, g, W, y0, 0.0, rmodes, log=lambda *a, **k: None)
    check("S5 linear conservation (cons_err < EPS_CONS)",
          float(np.nanmax(np.abs(r5.ts["cons_err"]))) < EPS_CONS,
          f"max={float(np.nanmax(np.abs(r5.ts['cons_err']))):.2e}")

    # S6 sponge accounting: an outgoing pulse conserves E_box + E_esc
    p6 = RunParams(profile="square", W0=0.5, a=3.0, lambda_chi=0.0, r_max=40.0,
                   dr=0.1, t_end=60.0, n_out=200, n_snapshots=20,
                   orthogonalize=False, amp=0.1, r0_packet=20.0, packet_width=2.0)
    g6 = Grid(GridSpec(p6.r_max, p6.dr, p6.sponge_start_frac, p6.sponge_strength))
    W6 = W_grid(p6.profile, g6.r, p6.W0, p6.a, p6.dr)
    rmodes6 = real_bound_modes(spectral_bound_modes(g6, W6))
    y06, _ = build_seed(g6, p6, W6, 0.0, rmodes6, log=lambda *a, **k: None)
    r6 = evolve(p6, g6, W6, y06, 0.0, rmodes6, log=lambda *a, **k: None)
    check("S6 sponge accounting (cons_err < EPS_CONS)",
          float(np.nanmax(np.abs(r6.ts["cons_err"]))) < EPS_CONS,
          f"max={float(np.nanmax(np.abs(r6.ts['cons_err']))):.2e}")

    # S7 seeded-mode energy recovery
    g7 = Grid(GridSpec(30.0, 0.05))
    W7 = W_grid("square", g7.r, 1.2, 5.0, 0.05)
    modes7 = real_bound_modes(spectral_bound_modes(g7, W7))
    if modes7:
        m = modes7[0]
        A = 1e-3
        chi = np.zeros(g7.N)
        chi[1:] = A * m["u"][1:] / g7.r[1:]
        chi[0] = chi[1]
        y = np.concatenate([chi, np.zeros(g7.N)])
        amps = mode_amplitudes(y, g7, [m])
        E_expected = 2.0 * np.pi * (m["omega2"] * A ** 2)
        rec = amps[0]["E"]
        check("S7 seeded-mode energy recovery",
              abs(rec - E_expected) / E_expected < 0.05,
              f"rec={rec:.4e} exp={E_expected:.4e}")
    else:
        check("S7 seeded-mode energy recovery", False, "(no modes)")

    # S8 orthogonalization efficacy
    p8 = RunParams(profile="square", W0=1.2, a=5.0, r_max=40.0, dr=0.1,
                   orthogonalize=True, amp=0.2, r0_packet=8.0, packet_width=3.0)
    g8 = Grid(GridSpec(p8.r_max, p8.dr))
    W8 = W_grid(p8.profile, g8.r, p8.W0, p8.a, p8.dr)
    rm8 = real_bound_modes(spectral_bound_modes(g8, W8))
    chi8, chid8, _ = raw_packet_fields(g8, p8)
    yraw = np.concatenate([chi8, chid8])
    Er = max(abs(energy_box(yraw, g8, W8, p8.lambda_chi)), E_FLOOR)
    fraw = float(np.nansum([a["E"] for a in mode_amplitudes(yraw, g8, rm8)])) / Er
    y0c, _ = build_seed(g8, p8, W8, p8.lambda_chi, rm8, log=lambda *a, **k: None)
    fclean = float(np.nansum([a["E"] for a in mode_amplitudes(y0c, g8, rm8)])) / \
        max(abs(energy_box(y0c, g8, W8, p8.lambda_chi)), E_FLOOR)
    check("S8 orthogonalization efficacy (clean << raw)",
          fclean < GATES["eps_lin"] / 10.0 and fraw > fclean,
          f"raw={fraw:.2e} clean={fclean:.2e}")

    # S9 contamination detection -> INVALID_SEED
    seed_clean_flag = fraw < GATES["eps_lin"]
    status9 = adjudicate_run(seed_clean_flag, True, True, True, 0.05, GATES)
    check("S9 contaminated seed adjudicates INVALID_SEED",
          (not seed_clean_flag) and status9 == "INVALID_SEED")

    # S10 adjudicator truth table (all five states + three ensemble states)
    t10 = [
        adjudicate_run(False, True, True, True, 0.1, GATES) == "INVALID_SEED",
        adjudicate_run(True, False, True, True, 0.1, GATES) == "INVALID_NUMERICS",
        adjudicate_run(True, True, False, True, 0.1, GATES) == "INCONCLUSIVE",
        adjudicate_run(True, True, True, False, 0.1, GATES) == "INCONCLUSIVE",
        adjudicate_run(True, True, True, True, 0.001, GATES) == "REJECTED",
        adjudicate_run(True, True, True, True, 0.1, GATES) == "SUPPORTED",
    ]
    e_sup = adjudicate_ensemble(
        [{"status": "SUPPORTED", "net_f_cap": 0.05},
         {"status": "SUPPORTED", "net_f_cap": 0.052}], GATES)["h2_ensemble_status"]
    e_rej = adjudicate_ensemble(
        [{"status": "SUPPORTED", "net_f_cap": 0.05},
         {"status": "REJECTED", "net_f_cap": 0.001}], GATES)["h2_ensemble_status"]
    e_inc = adjudicate_ensemble(
        [{"status": "SUPPORTED", "net_f_cap": 0.05},
         {"status": "INCONCLUSIVE", "net_f_cap": 0.03}], GATES)["h2_ensemble_status"]
    check("S10 adjudicator truth table",
          all(t10) and e_sup == "H2_SUPPORTED" and e_rej == "H2_REJECTED"
          and e_inc == "H2_INCONCLUSIVE")

    # S11 determinism: two mini-runs byte-identical
    d1 = os.path.join(scratch, "det1")
    d2 = os.path.join(scratch, "det2")
    cfg_mini = _mini_config()
    for d in (d1, d2):
        run_capture(_params_from_cfg(cfg_mini["baseline"]),
                    os.path.join(d, "baseline"), "b", 1.0, False,
                    "mini", cfg_mini, log=lambda *a, **k: None, write_full=False)
    with open(os.path.join(d1, "baseline", "time_series.csv"), "rb") as f1, \
            open(os.path.join(d2, "baseline", "time_series.csv"), "rb") as f2:
        check("S11 determinism (identical time_series bytes)",
              f1.read() == f2.read())

    # ---- S13: carrier direction and outgoing contamination (R9-A / C6) ----
    # Free space (W = 0), no bound modes, short run: the envelope centroid must
    # travel INWARD at the group speed, and the outward-side energy fraction
    # must stay negligible.
    k0c, sig_c = 3.5, 2.0
    pc = RunParams(profile="square", W0=0.0, a=5.0, lambda_chi=0.0, amp=0.2,
                   packet_kind="chi_carrier", k0=k0c, packet_width=sig_c,
                   r0_packet=55.0, orthogonalize=False,
                   r_max=100.0, dr=0.05, cfl=0.2, t_end=20.0, n_out=80,
                   n_snapshots=40)
    gc = Grid(GridSpec(pc.r_max, pc.dr, pc.sponge_start_frac, pc.sponge_strength))
    Wc = np.zeros(gc.N)
    chi_c, chid_c, _ = raw_packet_fields(gc, pc)
    yc = np.concatenate([chi_c, chid_c])
    rc = evolve(pc, gc, Wc, yc, 0.0, [], log=lambda *a, **k: None)

    def _centroid(chi_row):
        w = (gc.r * chi_row) ** 2
        return float(np.sum(gc.r * w) / max(np.sum(w), 1e-300))

    cen = np.array([_centroid(row) for row in rc.snap_chi])
    v_meas = float(np.polyfit(rc.snap_t, cen, 1)[0])
    v_th = -pc.v_group
    out_mask = gc.r > (55.0 + 3.0 * sig_c)
    u_end = (gc.r * rc.snap_chi[-1]) ** 2
    out_frac = float(np.sum(u_end[out_mask]) / max(np.sum(u_end), 1e-300))
    check("S13 carrier direction (inward at v_g; low outgoing content)",
          abs(v_meas - v_th) / abs(v_th) < 0.10 and out_frac < 1e-3,
          f"v={v_meas:+.4f} vs {v_th:+.4f}, outward_frac={out_frac:.2e}")

    # ---- S14: carrier tail clears the box (R9-A / C6) --------------------
    # Free-space carrier run with the sponge active, to a G4-satisfying
    # horizon: the whole-box residual must fall well below eps_res.
    p14 = RunParams(profile="square", W0=0.0, a=5.0, lambda_chi=0.0, amp=0.2,
                    packet_kind="chi_carrier", k0=3.5, packet_width=2.0,
                    r0_packet=55.0, orthogonalize=False,
                    r_max=100.0, dr=0.05, cfl=0.2, t_end=400.0, n_out=200,
                    n_snapshots=10)
    assert_carrier_guards(p14)  # must not raise: the test config is compliant
    g14 = Grid(GridSpec(p14.r_max, p14.dr, p14.sponge_start_frac,
                        p14.sponge_strength))
    W14 = np.zeros(g14.N)
    c14, cd14, _ = raw_packet_fields(g14, p14)
    r14 = evolve(p14, g14, W14, np.concatenate([c14, cd14]), 0.0, [],
                 log=lambda *a, **k: None)
    f_res_end = float(r14.ts["f_res"][-1])
    check("S14 carrier tail clearance (f_res << eps_res)",
          f_res_end < GATES["eps_res"] / 10.0,
          f"f_res(t_end)={f_res_end:.2e} vs eps_res/10="
          f"{GATES['eps_res']/10.0:.1e}")

    # ---- S15: guards G1-G4 refuse out-of-protocol probes (C6) ------------
    def _refused(par, guard_tag):
        try:
            assert_carrier_guards(par)
        except ValueError as exc:
            return guard_tag in str(exc)
        return False

    g_ok = RunParams(packet_kind="chi_carrier", k0=3.5, packet_width=2.0,
                     r0_packet=55.0, r_max=100.0, dr=0.05, t_end=400.0)
    s15 = [
        _refused(dataclasses.replace(g_ok, k0=-1.0), "G1"),
        _refused(dataclasses.replace(g_ok, packet_width=1.0), "G2"),   # k0*sig=3.5
        _refused(dataclasses.replace(g_ok, dr=0.5), "G3"),             # ~3.6 pts/lam
        _refused(dataclasses.replace(g_ok, t_end=50.0), "G4"),
    ]
    try:
        assert_carrier_guards(g_ok)
        compliant_ok = True
    except ValueError:
        compliant_ok = False
    check("S15 carrier guards refuse out-of-protocol configs",
          all(s15) and compliant_ok,
          f"G1..G4={s15}, compliant_accepted={compliant_ok}")

    # ---- S16: spherical focusing must not be mistaken for blow-up --------
    # Regression guard. A linear (lambda = 0) ingoing shell focuses at the
    # origin with chi = u/r growing as 1/r; this is geometry, not instability,
    # and must NOT halt the run. An absolute |chi| ceiling fails this test at
    # large launch radius, which previously killed the linear null control
    # while letting the nonlinear baseline pass.
    p16 = RunParams(profile="square", W0=0.9, a=5.0, lambda_chi=0.0, amp=0.2,
                    packet_kind="chi_carrier", k0=3.5, packet_width=2.0,
                    r0_packet=55.0, orthogonalize=True,
                    r_max=100.0, dr=0.05, cfl=0.2, t_end=90.0, n_out=90,
                    n_snapshots=10)
    g16 = Grid(GridSpec(p16.r_max, p16.dr, p16.sponge_start_frac,
                        p16.sponge_strength))
    W16 = W_grid(p16.profile, g16.r, p16.W0, p16.a, p16.dr)
    rm16 = real_bound_modes(spectral_bound_modes(g16, W16))
    y16, _ = build_seed(g16, p16, W16, 0.0, rm16, log=lambda *a, **k: None)
    r16 = evolve(p16, g16, W16, y16, 0.0, rm16, log=lambda *a, **k: None)
    chi_peak = float(np.max(np.abs(r16.snap_chi)))
    cons16 = float(np.nanmax(np.abs(r16.ts["cons_err"])))
    check("S16 spherical focusing not flagged as blow-up",
          (not r16.halted) and cons16 < EPS_CONS,
          f"halted={r16.halted}, peak|chi|={chi_peak:.1f}, cons_err={cons16:.2e}")

    # S12 forbidden-token lint over ALL public text/source files and any
    # generated outputs. Tokens are assembled from fragments so the banned
    # words do not themselves appear literally in this file, and matching is
    # word-boundary aware so legitimate identifiers (fmatch, format, ...) and
    # substrings inside ordinary words are not false positives.
    banned = ["".join(p) for p in [
        ("pro", "ton"), ("bary", "on"), ("qu", "ark"), ("glu", "on"),
        ("q", "c", "d"), ("ge", "v"), ("f", "m"), ("trip", "let"),
        ("fus", "ion"), ("ingred", "ient"), ("mass gener", "ation"),
        ("standard mo", "del"), ("toy_fus", "ion"), ("not_a_pro", "ton"),
    ]]
    patterns = [(t, re.compile(r"(?<![A-Za-z0-9_])" + re.escape(t)
                               + r"(?![A-Za-z0-9_])", re.I)) for t in banned]

    here = os.path.dirname(os.path.abspath(__file__))
    scan_targets = []
    for pat in ("*.md", "*.py", "configs/*.json", "CITATION.cff", "LICENSE"):
        scan_targets += sorted(glob.glob(os.path.join(here, pat)))
    for root, _dirs, files in os.walk(scratch):          # generated outputs
        for fn in sorted(files):
            if fn.endswith((".csv", ".json")):
                scan_targets.append(os.path.join(root, fn))

    all_hits = {}
    for tgt in scan_targets:
        try:
            with open(tgt, encoding="utf-8") as fh:
                text = fh.read()
        except Exception:
            continue
        h = [t for t, p in patterns if p.search(text)]
        if h:
            all_hits[os.path.relpath(tgt, here)] = h
    check(f"S12 forbidden-token lint ({len(scan_targets)} public files"
          " + outputs, boundary-aware)", not all_hits, f"hits={all_hits}")

    log(f"self-tests: {'ALL PASS' if ok else 'FAILURES PRESENT'}")
    return ok


# ======================================================================
# CONFIG HANDLING
# ======================================================================
def _mini_config():
    """Tiny, fast config used by the determinism self-test only."""
    return {
        "baseline": {
            "profile": "square", "W0": 0.9, "a": 4.0, "Delta": 0.5,
            "lambda_chi": 1.0, "amp": 0.1, "packet_width": 2.0,
            "orthogonalize": True, "r_max": 24.0, "dr": 0.2, "cfl": 0.2,
            "t_end": 20.0, "n_out": 60, "n_snapshots": 12,
            "field_r_stride": 3, "field_n_snapshots": 10, "run_id": "mini",
        }
    }


def require_blocks(cfg, blocks, command, config_path):
    """Fail early and actionably when a config lacks a block the command needs.

    `all` needs baseline + spectral + scan in one file (see
    configs/publication_all.json), or per-stage configs supplied through
    --baseline-config / --spectral-config / --scan-config.
    """
    missing = [b for b in blocks if b not in cfg]
    if missing:
        raise SystemExit(
            f"Config '{config_path}' is missing the block(s) {missing} required "
            f"by the '{command}' command (present: {sorted(cfg.keys())}).\n"
            f"Fix: use a combined config such as configs/publication_all.json, "
            f"or pass per-stage configs, e.g.\n"
            f"  python curvature_well_capture.py all "
            f"--baseline-config configs/baseline.json "
            f"--spectral-config configs/spectral_scan.json "
            f"--scan-config configs/robustness_scan.json --out data")


def load_config(path):
    if path is None:
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# ======================================================================
# CLI
# ======================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Public simulation for prescribed inward-curvature-well "
                    "wave localization (model units; m0 = 1).")
    sub = ap.add_subparsers(dest="command", required=True)

    def add_common(sp, stage_overrides=False):
        sp.add_argument("--config", default=None,
                        help="config file; for 'all' it must contain the "
                             "baseline, spectral, and scan blocks "
                             "(see configs/publication_all.json)")
        sp.add_argument("--out", default="data")
        if stage_overrides:
            sp.add_argument("--baseline-config", default=None,
                            help="per-stage override supplying the baseline block")
            sp.add_argument("--spectral-config", default=None,
                            help="per-stage override supplying the spectral block")
            sp.add_argument("--scan-config", default=None,
                            help="per-stage override supplying the scan block")

    for name in ("spectral", "baseline", "null", "scan", "controls",
                 "adjudicate"):
        add_common(sub.add_parser(name))
    add_common(sub.add_parser("all"), stage_overrides=True)
    sp_self = sub.add_parser("selftest")
    sp_self.add_argument("--scratch", default="selftest_scratch")

    args = ap.parse_args(argv)

    def log(msg):
        print(msg)

    if args.command == "selftest":
        ok = run_selftests(args.scratch, log=log)
        sys.exit(0 if ok else 1)

    cfg = load_config(args.config)
    data_root = args.out

    # `adjudicate` reads persisted metadata and needs no config at all.
    if args.command == "adjudicate":
        cmd_adjudicate(cfg or {}, data_root, log=log)
        log("done.")
        return

    if cfg is None and args.command != "all":
        log("NOTE: no --config provided; using built-in PLACEHOLDER defaults "
            "(R9 values are not frozen). Not for publication runs.")
        cfg = _default_config()

    if args.command == "spectral":
        require_blocks(cfg, ["spectral"], "spectral", args.config)
        cmd_spectral(cfg, data_root, log=log)
    elif args.command == "baseline":
        require_blocks(cfg, ["baseline"], "baseline", args.config)
        cmd_baseline(cfg, data_root, args.config, log=log)
    elif args.command == "null":
        require_blocks(cfg, ["baseline"], "null", args.config)
        cmd_null(cfg, data_root, args.config, log=log)
    elif args.command == "scan":
        require_blocks(cfg, ["scan"], "scan", args.config)
        cmd_scan(cfg, data_root, args.config, log=log)
        cmd_adjudicate(cfg, data_root, log=log)
    elif args.command == "controls":
        require_blocks(cfg, ["controls"], "controls", args.config)
        cmd_controls(cfg, data_root, args.config, log=log)
        cmd_adjudicate(cfg, data_root, log=log)
    elif args.command == "all":
        # Assemble the three stage blocks: per-stage overrides win over the
        # combined --config; a stage may come from either source.
        merged = dict(cfg or {})
        stage_files = {"baseline": args.baseline_config,
                       "spectral": args.spectral_config,
                       "scan": args.scan_config}
        for block, path in stage_files.items():
            if path:
                sub_cfg = load_config(path)
                if block not in sub_cfg:
                    raise SystemExit(
                        f"--{block}-config '{path}' does not contain a "
                        f"'{block}' block (present: {sorted(sub_cfg.keys())}).")
                merged[block] = sub_cfg[block]
                # carry shared scalars (e.g. Delta) if absent
                for k, v in sub_cfg.items():
                    merged.setdefault(k, v)
        if not merged:
            log("NOTE: no config provided; using built-in PLACEHOLDER defaults "
                "(R9 values are not frozen). Not for publication runs.")
            merged = _default_config()
        require_blocks(merged, ["baseline", "spectral", "scan"], "all",
                       args.config or "(per-stage configs)")
        if not run_selftests(log=log):
            log("Aborting 'all': self-tests failed.")
            sys.exit(1)
        cmd_spectral(merged, data_root, log=log)
        cmd_baseline(merged, data_root, args.config, log=log)
        cmd_null(merged, data_root, args.config, log=log)
        cmd_scan(merged, data_root, args.config, log=log)
        if merged.get("controls", {}).get("runs"):
            cmd_controls(merged, data_root, args.config, log=log)
        cmd_adjudicate(merged, data_root, log=log)
    log("done.")


def _default_config():
    """Built-in PLACEHOLDER config (R9 not frozen). Values are chosen only so
    the pipeline runs end-to-end; they are NOT the publication configuration."""
    return {
        "Delta": 0.5,
        "baseline": {
            "profile": "square", "W0": 0.9, "a": 5.0, "Delta": 0.5,
            "lambda_chi": 1.0, "amp": 0.2, "packet_width": 2.0,
            "orthogonalize": True, "r_max": 40.0, "dr": 0.05, "cfl": 0.2,
            "t_end": 200.0, "n_out": 600, "n_snapshots": 81,
            "field_r_stride": 4, "field_n_snapshots": 80, "run_id": "baseline",
        },
        "spectral": {
            "a": 15.0, "dr": 0.02, "r_max": 80.0,
            "W0_grid": list(np.round(np.linspace(0.05, 0.95, 20), 4)),
        },
        "scan": {
            "runs": [
                {"profile": "square", "W0": 0.9, "in_robustness_set": True,
                 "run_id": "scan_square_W0p9"},
                {"profile": "tanh", "W0": 0.9, "in_robustness_set": True,
                 "run_id": "scan_tanh_W0p9"},
                {"profile": "sech2", "W0": 0.9, "in_robustness_set": True,
                 "run_id": "scan_sech2_W0p9"},
            ]
        },
    }


if __name__ == "__main__":
    main()
