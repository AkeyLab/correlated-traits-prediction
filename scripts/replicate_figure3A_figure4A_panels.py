#!/usr/bin/env python3
"""Independent replications of Figure 3 panel A and Figure 4 panel A, with error bars.

Same purpose as ``replicate_figure3A_column.py`` but over the full grid: 21
environmental correlations by 9 helper heritabilities. Per-replicate values are
retained so every cell gets a Monte Carlo standard error.

Those standard errors are what makes Figure 2D-F possible. The published runs
stored only the 30-replicate averages, so the uncertainty plotted in Figure 2D-F
comes from these replications -- identical design (same sample size, same number
of replicates, same causal-variant structure), different seed.

    Figure 3A design   h2_t = 0.9   ->  Figure 2D and 2E
    Figure 4A design   h2_t = 0.1   ->  Figure 2F

Both designs fix rho_g = 0.9 and alpha = 0.1. Each design starts its own
generator from REPLICATION_SEED, so the Figure 3A replication is unchanged by
the addition of the Figure 4A one.

The explicit least-squares route is not repeated here. It was already shown in
``replicate_figure3A_column.py`` to agree with the closed-form squared multiple
correlation to about 1e-14, so only the closed form is evaluated.

Runtime: roughly 10-30 minutes per design on one core.
Output:  tables/replication_fig3A_panel.npz
         tables/replication_fig4A_panel.npz
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

RHO_G = 0.9
ALPHA = 0.1

#: (output file, target heritability)
DESIGNS = [
    ("replication_fig3A_panel.npz", 0.9),
    ("replication_fig4A_panel.npz", 0.1),
]


def replicate(config, h2_target):
    """Per-replicate gains over the panel-A grid, shape (reps, rho_e, h2_h)."""
    rng = np.random.default_rng(REPLICATION_SEED)
    gains = np.zeros((config.n_replicates, len(RHO_GRID), len(H2_HELPER_GRID)))

    for replicate_index in range(config.n_replicates):
        genotypes = draw_genotypes(rng, config)
        a, b, c = draw_basis(rng, config)
        g_a, g_b, g_c = genotypes @ a, genotypes @ b, genotypes @ c
        for i, rho_e in enumerate(RHO_GRID):
            for j, h2_helper in enumerate(H2_HELPER_GRID):
                gains[replicate_index, i, j] = relative_gain(
                    rng, config, g_a, g_b, g_c,
                    h2_target, h2_helper, RHO_G, rho_e, ALPHA,
                )
        print("  replicate %d/%d done" % (replicate_index + 1, config.n_replicates),
              flush=True)
    return gains


def main() -> None:
    paths.ensure_output_dirs()
    config = SimulationConfig()

    for filename, h2_target in DESIGNS:
        print("=== h2_t = %.1f -> %s ===" % (h2_target, filename))
        gains = replicate(config, h2_target)

        mean = gains.mean(axis=0)
        standard_deviation = gains.std(axis=0, ddof=1)
        standard_error = standard_deviation / np.sqrt(config.n_replicates)

        output = paths.TABLES / filename
        np.savez(
            output,
            mean=mean, sd=standard_deviation, se=standard_error,
            rho_e_grid=RHO_GRID, h2_h_grid=H2_HELPER_GRID,
            reps=config.n_replicates, n=config.n_individuals,
            h2_t=h2_target, rho_g=RHO_G, alpha=ALPHA,
        )
        print("Wrote %s" % output)
        print("Replication mean gain range: %.4f to %.4f" % (mean.min(), mean.max()))
        print("Per-cell standard error range: %.5f to %.5f\n"
              % (standard_error.min(), standard_error.max()))


if __name__ == "__main__":
    main()
