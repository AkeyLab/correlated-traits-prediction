#!/usr/bin/env python3
"""Independent replication of one column of Figure 3 panel A.

The published heatmaps average 30 replicates per cell but report no
uncertainty, so a reader cannot tell whether a small gap between simulation and
the closed form is a real discrepancy or just Monte Carlo noise. This script
re-runs the same design with a *different seed*, keeping per-replicate values,
and answers two questions the published figure cannot:

1.  **Are the deviations consistent with Monte Carlo error?**
    Retaining per-replicate values gives a standard error per cell and hence a
    z-score against the closed-form prediction.

2.  **Does an explicit regression agree with the closed-form identity?**
    The Methods describe fitting joint models. The figure-generating code
    instead evaluates the algebraically equivalent squared multiple
    correlation. This script computes both and reports the largest discrepancy;
    it comes out around 1e-14, i.e. floating-point noise.

The column analysed is ``h2_h = 0.6`` at the Figure 3 parameters
(``h2_t = 0.9``, ``rho_g = 0.9``, ``alpha = 0.1``), sweeping the environmental
correlation from -1 to 1.

This is a replication, not a reproduction: the seed differs from the published
run on purpose.

Runtime: a few minutes.
Output:  tables/replication_fig3A_column_hh2_0.6.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import REPLICATION_SEED, SimulationConfig, paths  # noqa: E402
from correlated_traits.simulate import (  # noqa: E402
    draw_basis, draw_genotypes, pearson, simulate_cell,
)
from correlated_traits.theory import gain_closed_form, joint_r2  # noqa: E402

H2_TARGET = 0.9
RHO_G = 0.9
ALPHA = 0.1
R2_BASE_NOMINAL = ALPHA * H2_TARGET
COLUMN_H2_HELPER = 0.6
RHO_E_GRID = np.linspace(-1, 1, 21)


def ordinary_least_squares_r2(response, predictors):
    """In-sample R^2 from a least-squares fit with an intercept."""
    design = np.column_stack([np.ones(len(response))] + list(predictors))
    coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ coefficients
    centred = response - response.mean()
    return 1.0 - (residual @ residual) / (centred @ centred)


def main() -> None:
    paths.ensure_output_dirs()
    config = SimulationConfig()
    rng = np.random.default_rng(REPLICATION_SEED)

    n_rho = len(RHO_E_GRID)
    gains_closed_form = np.zeros((config.n_replicates, n_rho))
    gains_least_squares = np.zeros((config.n_replicates, n_rho))
    r2_base_empirical = np.zeros((config.n_replicates, n_rho))

    for replicate in range(config.n_replicates):
        genotypes = draw_genotypes(rng, config)
        a, b, c = draw_basis(rng, config)
        g_a, g_b, g_c = genotypes @ a, genotypes @ b, genotypes @ c

        for i, rho_e in enumerate(RHO_E_GRID):
            target, helper, baseline = simulate_cell(
                rng, config, g_a, g_b, g_c,
                H2_TARGET, COLUMN_H2_HELPER, RHO_G, rho_e, ALPHA,
            )

            r_tb = pearson(target, baseline)
            r_th = pearson(target, helper)
            r_hb = pearson(helper, baseline)
            r2_base = r_tb ** 2
            gains_closed_form[replicate, i] = (joint_r2(r_tb, r_th, r_hb) - r2_base) / r2_base

            r2_base_fit = ordinary_least_squares_r2(target, [baseline])
            r2_joint_fit = ordinary_least_squares_r2(target, [baseline, helper])
            gains_least_squares[replicate, i] = (r2_joint_fit - r2_base_fit) / r2_base_fit

            r2_base_empirical[replicate, i] = r2_base

    predicted = gain_closed_form(
        H2_TARGET, COLUMN_H2_HELPER, RHO_G, RHO_E_GRID, R2_BASE_NOMINAL
    )
    observed = gains_closed_form.mean(axis=0)
    standard_error = gains_closed_form.std(axis=0, ddof=1) / np.sqrt(config.n_replicates)
    z_scores = (observed - predicted) / standard_error

    print("Independent replication: %d replicates, n = %d per replicate"
          % (config.n_replicates, config.n_individuals))
    print("Column h2_h = %.1f; h2_t = %.1f, rho_g = %.1f, alpha = %.1f, "
          "nominal R2_base = %g"
          % (COLUMN_H2_HELPER, H2_TARGET, RHO_G, ALPHA, R2_BASE_NOMINAL))
    print()
    print("Closed form vs explicit least squares: max |difference| over all %d "
          "cell-replicates = %.3e"
          % (config.n_replicates * n_rho,
             np.max(np.abs(gains_closed_form - gains_least_squares))))
    print("Empirical R2_base: mean = %.5f (nominal %g), sd = %.5f"
          % (r2_base_empirical.mean(), R2_BASE_NOMINAL, r2_base_empirical.std(ddof=1)))
    print()
    print("%7s %9s %13s %8s %9s %7s"
          % ("rho_e", "theory", "replic. mean", "SE", "diff", "z"))
    for i, rho_e in enumerate(RHO_E_GRID):
        print("%7.1f %9.4f %13.4f %8.4f %9.4f %7.2f"
              % (rho_e, predicted[i], observed[i], standard_error[i],
                 observed[i] - predicted[i], z_scores[i]))

    print()
    print("|z| summary: max = %.2f, mean = %.2f, cells with |z| > 2: %d of %d"
          % (np.max(np.abs(z_scores)), np.mean(np.abs(z_scores)),
             int(np.sum(np.abs(z_scores) > 2)), n_rho))
    print("Mean signed relative deviation across the column = %+.3f%%"
          % (100 * np.mean((observed - predicted) / predicted)))
    print("Typical per-cell Monte Carlo relative SE = %.2f%%"
          % (100 * np.mean(standard_error / predicted)))

    output = paths.TABLES / "replication_fig3A_column_hh2_0.6.csv"
    with open(output, "w") as handle:
        handle.write("rho_e,theory_R2gain,replication_mean,replication_se,z\n")
        for i, rho_e in enumerate(RHO_E_GRID):
            handle.write("%.1f,%.6f,%.6f,%.6f,%.4f\n"
                         % (rho_e, predicted[i], observed[i],
                            standard_error[i], z_scores[i]))
    print("\nWrote %s" % output)


if __name__ == "__main__":
    main()
