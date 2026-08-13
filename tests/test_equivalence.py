"""Regression test: the package must be bit-for-bit identical to the original scripts.

The figures in the manuscript were produced by a set of standalone scripts in
which the simulation, the parameter grids and the plotting were interleaved in
one file. This repository factors that into ``correlated_traits``. Refactoring
simulation code is dangerous, because results depend not only on the arithmetic
but on the *order* in which random numbers are drawn -- a change that looks
harmless can silently move every number in a figure.

This test pins that down. It contains a verbatim transcription of the original
implementation and asserts that the package reproduces it exactly, using the
same seed at reduced scale so the test runs in seconds rather than an hour.

Reduced scale is safe for this purpose: an ordering mistake changes the very
first cell, so it does not hide behind a small sample size. The full-scale
check against the published arrays is a separate script,
``scripts/verify_reproducibility.py``.

Run with:  python -m pytest tests/ -v      (or simply: python tests/test_equivalence.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import SimulationConfig  # noqa: E402
from correlated_traits.simulate import heatmap_panel  # noqa: E402

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


def test_package_matches_original_bit_for_bit():
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


def test_closed_form_matches_simulation_identity():
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


def test_closed_form_matches_the_published_arrays():
    """The paper's headline equation must reproduce the simulated figures.

    ``gain_closed_form`` is Equation 12 of the manuscript and the expression the
    public calculator implements, yet nothing here exercised it. This checks it
    against the published simulation arrays: the closed form is a population
    quantity and the arrays are 30-replicate Monte Carlo averages, so they agree
    closely but not exactly, and the tolerance below is on the correlation across
    a whole panel rather than cell by cell.
    """
    import pathlib

    from correlated_traits.theory import gain_closed_form

    rho = np.linspace(-1, 1, 21)
    h2h = np.round(np.linspace(0.1, 0.9, 9), 2)
    grid_rho, grid_h2h = np.meshgrid(rho, h2h, indexing="ij")

    reference = pathlib.Path(__file__).resolve().parents[1] / "reference" / "arrays"
    if not reference.is_dir():
        raise AssertionError("reference/arrays is missing; cannot validate the closed form")

    # (file stem, h2_t, alpha, panel) -- panel A fixes rho_g = 0.9 and varies rho_e,
    # panel B fixes rho_e = 0.9 and varies rho_g.
    cases = [
        ("Fig3_scenario1_highh", 0.9, 0.1),
        ("Fig4_scenario2_lowh", 0.1, 0.1),
        ("S2_scenario3_lowh_goodbase", 0.1, 0.9),
        ("S3_scenario4_highh_goodbase", 0.9, 0.9),
    ]

    checked = 0
    for stem, h2_t, alpha in cases:
        r2_base = alpha * h2_t
        for panel in ("A", "B"):
            path = reference / ("%s_%s.npy" % (stem, panel))
            if not path.exists():
                raise AssertionError("missing reference array %s" % path)
            simulated = np.load(path)
            if panel == "A":
                predicted = gain_closed_form(h2_t, grid_h2h, 0.9, grid_rho, r2_base)
            else:
                predicted = gain_closed_form(h2_t, grid_h2h, grid_rho, 0.9, r2_base)

            r = np.corrcoef(predicted.ravel(), simulated.ravel())[0, 1]
            assert r > 0.99, "%s panel %s: closed form vs simulation r = %.4f" % (stem, panel, r)
            checked += 1

    assert checked == 8, "expected to check 8 panels, checked %d" % checked


def test_phenotypic_correlation_stays_a_correlation():
    """rho_p must lie in [-1, 1] for every admissible parameter combination."""
    from correlated_traits.theory import phenotypic_correlation

    rng = np.random.default_rng(7)
    h2_t = rng.uniform(0.01, 0.99, 20000)
    h2_h = rng.uniform(0.01, 0.99, 20000)
    rho_g = rng.uniform(-1, 1, 20000)
    rho_e = rng.uniform(-1, 1, 20000)
    rho_p = phenotypic_correlation(h2_t, h2_h, rho_g, rho_e)
    assert np.all(np.abs(rho_p) <= 1 + 1e-12), "rho_p left [-1,1]: max |rho_p| = %.6f" % np.abs(rho_p).max()


if __name__ == "__main__":
    test_package_matches_original_bit_for_bit()
    print("PASS  package reproduces the original implementation bit-for-bit")
    test_closed_form_matches_simulation_identity()
    print("PASS  two-predictor R^2 closed form matches an explicit least-squares fit")
    test_closed_form_matches_the_published_arrays()
    print("PASS  gain_closed_form reproduces all eight published simulation panels")
    test_phenotypic_correlation_stays_a_correlation()
    print("PASS  phenotypic correlation stays within [-1, 1]")
