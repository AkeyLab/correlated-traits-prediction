#!/usr/bin/env python3
"""Figures S2 and S3 -- the same heatmaps under a STRONG baseline.

    Figure S3   Scenario 3: weakly  heritable target (h2_t = 0.1), alpha = 0.9
    Figure S2   Scenario 4: highly heritable target (h2_t = 0.9), alpha = 0.9

The only change from Figures 3 and 4 is ``alpha``: the baseline predictor now
captures 90% of the target heritability instead of 10%. These two figures are
what the manuscript uses to argue that helper traits still pay off for a weakly
heritable target even when the baseline is already strong, but stop paying off
for a highly heritable one.

Note the file naming, which is a known trap: the stems follow the SCENARIO
numbering, not the manuscript's figure numbers, and the two are off by one.
``S2_scenario3_*`` prints in the paper as Figure S3; ``S3_scenario4_*`` prints
as Figure S2. The manuscript numbers supplementary figures by order of
appearance in Additional file 1, where the Scenario 4 panel comes first.

Runtime is comparable to Figures 3 and 4, roughly 30-60 minutes on one core.

IMPORTANT -- as with Figures 3 and 4, both figures come from a single random
generator in a fixed order (S2 panel A, S2 panel B, S3 panel A, S3 panel B).
Do not split or reorder them.

Outputs (under results/):
    arrays/S2_scenario3_lowh_goodbase_{A,B}.npy
    arrays/S3_scenario4_highh_goodbase_{A,B}.npy
    figures/S2_scenario3_lowh_goodbase.png
    figures/S3_scenario4_highh_goodbase.png
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

ALPHA_STRONG_BASELINE = 0.9

FIGURES = [
    ("S2_scenario3_lowh_goodbase", 0.1, "Figure S2  Scenario 3: weakly heritable target"),
    ("S3_scenario4_highh_goodbase", 0.9, "Figure S3  Scenario 4: highly heritable target"),
]


def main() -> None:
    paths.ensure_output_dirs()
    config = SimulationConfig()
    rng = np.random.default_rng(HEATMAP_SEED)

    for stem, h2_target, label in FIGURES:
        print("=== %s (h2_t = %.1f, alpha = %.1f) ===" % (label, h2_target, ALPHA_STRONG_BASELINE))
        panels = {}
        for mode in ("A", "B"):
            panels[mode] = heatmap_panel(
                rng, config, h2_target, ALPHA_STRONG_BASELINE, mode, RHO_GRID, H2_HELPER_GRID
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
