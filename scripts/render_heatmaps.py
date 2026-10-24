#!/usr/bin/env python3
"""Redraw Figures 3, 4, S1 and S2 from saved arrays, without re-simulating.

Reads results/arrays/ if present, otherwise the published reference/arrays/,
so the house style can be applied in seconds. The values are identical to the
simulation scripts' output (scripts/verify_reproducibility.py checks that).

Outputs (under results/figures/):
    Fig3_scenario1_highh.{png,pdf}         Figure 3
    Fig4_scenario2_lowh.{png,pdf}          Figure 4
    S3_scenario4_highh_goodbase.{png,pdf}  Supplementary Figure 1
    S2_scenario3_lowh_goodbase.{png,pdf}   Supplementary Figure 2
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import H2_HELPER_GRID, RHO_GRID, paths  # noqa: E402
from correlated_traits.plotting import draw_heatmap_pair  # noqa: E402

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "arrays"
STEMS = ["Fig3_scenario1_highh", "Fig4_scenario2_lowh",
         "S3_scenario4_highh_goodbase", "S2_scenario3_lowh_goodbase"]


def load(name):
    for folder in (paths.ARRAYS, REFERENCE):
        if (folder / name).exists():
            return np.load(folder / name)
    raise SystemExit("missing array %s" % name)


def main() -> None:
    paths.ensure_output_dirs()
    for stem in STEMS:
        out = draw_heatmap_pair(load(stem + "_A.npy"), load(stem + "_B.npy"),
                                RHO_GRID, H2_HELPER_GRID, paths.FIGURES / stem)
        print("figure -> %s" % Path(out).name)


if __name__ == "__main__":
    main()
