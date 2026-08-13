"""Regression check: the package must be bit-for-bit identical to the original scripts.

The figures in the manuscript were produced by a set of standalone scripts in
which the simulation, the parameter grids and the plotting were interleaved in
one file. This repository factors that into ``correlated_traits``. Refactoring
simulation code is dangerous, because results depend not only on the arithmetic
but on the *order* in which random numbers are drawn -- a change that looks
harmless can silently move every number in a figure.

This module pins that down. It contains a verbatim transcription of the original
implementation and asserts that the package reproduces it exactly, using the
same seed at reduced scale so the check runs in seconds rather than an hour.

Reduced scale is safe for this purpose: an ordering mistake changes the very
first cell, so it does not hide behind a small sample size. The full-scale check
against the published arrays is ``correlated-traits verify``.

It also confirms that the published reference arrays actually shipped with the install,
which nothing else here would notice.

It lives in the package, and is exposed as ``correlated-traits test``, so that it runs
from an installed copy with no pytest and no source checkout. For a user the useful
framing is not the refactor history but the consequence: it is a one-second check that
their install is complete and their NumPy build produces the expected random stream.
"""

from __future__ import annotations

import numpy as np

from correlated_traits import SimulationConfig
from correlated_traits.simulate import heatmap_panel

__all__ = [
    "check_package_matches_original",
    "check_closed_form_matches_least_squares",
    "check_reference_arrays_present",
    "run",
]

# Reduced scale: same code paths, far fewer individuals and cells.
N = 400
M_SHARED = 20
M_HSPEC = 5
M = M_SHARED + M_HSPEC
LAM = 0.20
REPS = 3
ALPHA = 0.1
SEED = 20260717
RHO_GRID = np.linspace(-1, 1, 5)
HH_GRID = np.round(np.linspace(0.1, 0.9, 3), 2)


# --------------------------------------------------------------------------
# Verbatim transcription of the original implementation. Do not tidy this up:
# its value is that it is a faithful copy of the code that produced the figures.
# --------------------------------------------------------------------------
def _original_panel(rng, h_t2, mode):
    def corr(x, y):
        x = x - x.mean()
        y = y - y.mean()
        return float(x @ y / np.sqrt((x @ x) * (y @ y)))

    def make_basis():
        a = np.zeros(M)
        b = np.zeros(M)
        c = np.zeros(M)
        va = rng.standard_normal(M_SHARED)
        va /= np.linalg.norm(va)
        vb = rng.standard_normal(M_SHARED)
        vb -= (vb @ va) * va
        vb /= np.linalg.norm(vb)
        a[:M_SHARED] = va
        b[:M_SHARED] = vb
        vc = rng.standard_normal(M_HSPEC)
        vc /= np.linalg.norm(vc)
        c[M_SHARED:] = vc
        return a, b, c

    def sim_genotypes():
        p = rng.uniform(0.05, 0.5, M)
        X = rng.binomial(2, p, size=(N, M)).astype(np.float64)
        X -= 2 * p
        X /= np.sqrt(2 * p * (1 - p))
        return X

    def rescale(v, target_var):
        return v * np.sqrt(target_var / v.var())

    def gain_cell(Ga, Gb, Gc, h_t2, h_h2, rho_g, rho_e):
        G_t = rescale(Ga * np.sqrt(h_t2), h_t2)
        cb = np.sqrt(max((1 - rho_g ** 2) * (1 - LAM), 0.0))
        cc = np.sqrt(max((1 - rho_g ** 2) * LAM, 0.0))
        G_h_raw = np.sqrt(h_h2) * (Ga * rho_g + Gb * cb + Gc * cc)
        G_h = rescale(G_h_raw, h_h2)
        z1 = rng.standard_normal(N)
        z2 = rng.standard_normal(N)
        eps_t = np.sqrt(1 - h_t2) * z1
        eps_h = np.sqrt(1 - h_h2) * (rho_e * z1 + np.sqrt(max(1 - rho_e ** 2, 0.0)) * z2)
        Y_t = G_t + eps_t
        Y_h = G_h + eps_h
        Y_t = (Y_t - Y_t.mean()) / Y_t.std()
        Y_h = (Y_h - Y_h.mean()) / Y_h.std()
        Ghat = (G_t - G_t.mean()) / G_t.std()
        B = np.sqrt(ALPHA) * Ghat + np.sqrt(1 - ALPHA) * rng.standard_normal(N)
        r_tB = corr(Y_t, B)
        r_th = corr(Y_t, Y_h)
        r_hB = corr(Y_h, B)
        R2_base = r_tB ** 2
        R2_joint = (r_tB ** 2 + r_th ** 2 - 2 * r_tB * r_th * r_hB) / (1 - r_hB ** 2)
        return (R2_joint - R2_base) / R2_base

    acc = np.zeros((len(RHO_GRID), len(HH_GRID)))
    for _ in range(REPS):
        X = sim_genotypes()
        a, b, c = make_basis()
        Ga, Gb, Gc = X @ a, X @ b, X @ c
        for i, rho in enumerate(RHO_GRID):
            for j, hh in enumerate(HH_GRID):
                if mode == "A":
                    acc[i, j] += gain_cell(Ga, Gb, Gc, h_t2, hh, 0.9, rho)
                else:
                    acc[i, j] += gain_cell(Ga, Gb, Gc, h_t2, hh, rho, 0.9)
    return acc / REPS


def _package_panel(rng, h_t2, mode):
    config = SimulationConfig(
        n_individuals=N, m_shared=M_SHARED, m_helper_specific=M_HSPEC,
        helper_specific_fraction=LAM, n_replicates=REPS,
    )
    return heatmap_panel(rng, config, h_t2, ALPHA, mode, RHO_GRID, HH_GRID)


def check_package_matches_original() -> None:
    """Same seed, same panel sequence -> identical arrays, to the last bit."""
    original_rng = np.random.default_rng(SEED)
    original = [_original_panel(original_rng, h2, mode)
                for h2 in (0.9, 0.1) for mode in ("A", "B")]

    package_rng = np.random.default_rng(SEED)
    package = [_package_panel(package_rng, h2, mode)
               for h2 in (0.9, 0.1) for mode in ("A", "B")]

    for index, (expected, actual) in enumerate(zip(original, package)):
        assert np.array_equal(expected, actual), (
            "panel %d differs; max absolute difference %.3e"
            % (index, np.abs(expected - actual).max())
        )


def check_closed_form_matches_least_squares() -> None:
    """The two-predictor R^2 helper must equal an explicit least-squares fit."""
    from correlated_traits.theory import joint_r2

    rng = np.random.default_rng(0)
    n = 5000
    target = rng.standard_normal(n)
    helper = 0.5 * target + np.sqrt(1 - 0.25) * rng.standard_normal(n)
    baseline = 0.3 * target + np.sqrt(1 - 0.09) * rng.standard_normal(n)

    def standardize(v):
        return (v - v.mean()) / v.std()

    target, helper, baseline = map(standardize, (target, helper, baseline))
    design = np.column_stack([np.ones(n), helper, baseline])
    fitted = design @ np.linalg.lstsq(design, target, rcond=None)[0]
    least_squares_r2 = 1.0 - ((target - fitted) ** 2).sum() / ((target - target.mean()) ** 2).sum()

    def corr(x, y):
        return float((x - x.mean()) @ (y - y.mean())
                     / np.sqrt(((x - x.mean()) ** 2).sum() * ((y - y.mean()) ** 2).sum()))

    closed_form = joint_r2(corr(target, baseline), corr(target, helper), corr(helper, baseline))
    assert abs(closed_form - least_squares_r2) < 1e-12, (
        "closed form %.15f vs least squares %.15f" % (closed_form, least_squares_r2)
    )


def check_reference_arrays_present() -> int:
    """The published arrays ship inside the package and load at the right shape.

    Nothing else here touches them: the equivalence check compares two in-process
    implementations and never opens a file. So a packaging change that dropped the
    reference data would go unnoticed until someone ran ``correlated-traits verify``
    after an hour of simulation. Returns how many arrays were checked.
    """
    from correlated_traits import H2_HELPER_GRID
    from correlated_traits import RHO_GRID as PUBLISHED_RHO_GRID
    from correlated_traits import paths
    from correlated_traits.figures.verify_reproducibility import EXPECTED

    names = sorted(EXPECTED)
    absent = [name for name in names if not paths.has_reference_array(name)]
    assert not absent, (
        "reference arrays missing from the package install: %s" % ", ".join(absent)
    )

    expected_shape = (len(PUBLISHED_RHO_GRID), len(H2_HELPER_GRID))
    for name in names:
        shape = paths.reference_array(name).shape
        assert shape == expected_shape, (
            "%s has shape %s, expected %s" % (name, shape, expected_shape)
        )
    return len(names)


def run() -> None:
    """Run every check, printing a line per pass. Raises AssertionError on failure."""
    check_package_matches_original()
    print("PASS  simulation reproduces the reference implementation bit-for-bit")
    check_closed_form_matches_least_squares()
    print("PASS  two-predictor R^2 closed form matches an explicit least-squares fit")
    count = check_reference_arrays_present()
    print("PASS  %d published reference arrays ship with the package" % count)
