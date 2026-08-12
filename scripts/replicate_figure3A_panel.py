#!/usr/bin/env python3
"""Independent replication of the whole of Figure 3 panel A, with error bars.

Same purpose as ``replicate_figure3A_column.py`` but over the full grid: 21
environmental correlations by 9 helper heritabilities. Per-replicate values are
retained so every cell gets a Monte Carlo standard error.

Those standard errors are what makes Figure S1 possible. The published run
stored only the 30-replicate averages, so the uncertainty plotted in Figure S1
comes from this replication -- identical design (same sample size, same number
of replicates, same causal-variant structure), different seed.

The explicit least-squares route is not repeated here. It was already shown in
``replicate_figure3A_column.py`` to agree with the closed-form squared multiple
correlation to about 1e-14, so only the closed form is evaluated.

Runtime: roughly 10-30 minutes on one core.
Output:  tables/replication_fig3A_panel.npz
         keys: mean, sd, se, rho_e_grid, h2_h_grid, and the design parameters
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import (  # noqa: E402
    H2_HELPER_GRID, REPLICATION_SEED, RHO_GRID, SimulationConfig, paths,
)
from correlated_traits.simulate import (  # noqa: E402
    draw_basis, draw_genotypes, relative_gain,
)

H2_TARGET = 0.9
RHO_G = 0.9
ALPHA = 0.1


def main() -> None:
    paths.ensure_output_dirs()
    config = SimulationConfig()
    rng = np.random.default_rng(REPLICATION_SEED)

    gains = np.zeros((config.n_replicates, len(RHO_GRID), len(H2_HELPER_GRID)))

    for replicate in range(config.n_replicates):
        genotypes = draw_genotypes(rng, config)
        a, b, c = draw_basis(rng, config)
        g_a, g_b, g_c = genotypes @ a, genotypes @ b, genotypes @ c
        for i, rho_e in enumerate(RHO_GRID):
            for j, h2_helper in enumerate(H2_HELPER_GRID):
                gains[replicate, i, j] = relative_gain(
                    rng, config, g_a, g_b, g_c,
                    H2_TARGET, h2_helper, RHO_G, rho_e, ALPHA,
                )
        print("  replicate %d/%d done" % (replicate + 1, config.n_replicates), flush=True)

    mean = gains.mean(axis=0)
    standard_deviation = gains.std(axis=0, ddof=1)
    standard_error = standard_deviation / np.sqrt(config.n_replicates)

    output = paths.TABLES / "replication_fig3A_panel.npz"
    np.savez(
        output,
        mean=mean, sd=standard_deviation, se=standard_error,
        rho_e_grid=RHO_GRID, h2_h_grid=H2_HELPER_GRID,
        reps=config.n_replicates, n=config.n_individuals,
        h2_t=H2_TARGET, rho_g=RHO_G, alpha=ALPHA,
    )
    print("\nWrote %s" % output)
    print("Replication mean gain range: %.4f to %.4f" % (mean.min(), mean.max()))
    print("Per-cell standard error range: %.5f to %.5f"
          % (standard_error.min(), standard_error.max()))


if __name__ == "__main__":
    main()
