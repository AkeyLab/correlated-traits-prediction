#!/usr/bin/env python3
"""Figures 3 and 4 -- helper-trait characteristics under a POOR baseline.

    Figure 3   Scenario 1: highly heritable target (h2_t = 0.9), alpha = 0.1
    Figure 4   Scenario 2: weakly  heritable target (h2_t = 0.1), alpha = 0.1

Each figure has two panels: panel A fixes the genetic correlation at 0.9 and
sweeps the environmental correlation; panel B fixes the environmental
correlation at 0.9 and sweeps the genetic correlation.

Runtime is roughly 30-60 minutes on one core: 4 panels x 30 replicates x 189
grid cells, each cell simulating 20,000 individuals.

IMPORTANT -- do not split this script. Both figures are generated from a single
random generator, in the order Figure 3 panel A, Figure 3 panel B, Figure 4
panel A, Figure 4 panel B. Running them separately, or reordering them, changes
which random numbers each panel receives and the published arrays will no
longer reproduce bit-for-bit.

Outputs (under results/):
    arrays/Fig3_scenario1_highh_{A,B}.npy
    arrays/Fig4_scenario2_lowh_{A,B}.npy
    figures/Fig3_scenario1_highh.png
    figures/Fig4_scenario2_lowh.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import (  # noqa: E402
    H2_HELPER_GRID, HEATMAP_SEED, RHO_GRID, SimulationConfig, paths,
)
from correlated_traits.plotting import draw_heatmap_pair  # noqa: E402
from correlated_traits.simulate import heatmap_panel  # noqa: E402

ALPHA_POOR_BASELINE = 0.1

FIGURES = [
    ("Fig3_scenario1_highh", 0.9, "Figure 3  Scenario 1: highly heritable target"),
    ("Fig4_scenario2_lowh", 0.1, "Figure 4  Scenario 2: weakly heritable target"),
]


def main() -> None:
    paths.ensure_output_dirs()
    config = SimulationConfig()
    rng = np.random.default_rng(HEATMAP_SEED)

    for stem, h2_target, label in FIGURES:
        print("=== %s (h2_t = %.1f, alpha = %.1f) ===" % (label, h2_target, ALPHA_POOR_BASELINE))
        panels = {}
        for mode in ("A", "B"):
            panels[mode] = heatmap_panel(
                rng, config, h2_target, ALPHA_POOR_BASELINE, mode, RHO_GRID, H2_HELPER_GRID
            )
            array_path = paths.ARRAYS / ("%s_%s.npy" % (stem, mode))
            np.save(array_path, panels[mode])
            print("    panel %s  peak R2_gain = %8.3f   ->  %s"
                  % (mode, panels[mode].max(), array_path.name))

        figure_path = draw_heatmap_pair(
            panels["A"], panels["B"], RHO_GRID, H2_HELPER_GRID,
            paths.FIGURES / ("%s.png" % stem),
        )
        print("    figure -> %s\n" % figure_path.name)


if __name__ == "__main__":
    main()
