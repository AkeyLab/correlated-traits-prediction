#!/usr/bin/env python3
"""Supplementary Figure 4: validation and illustration of the multi-helper gain.

    panel A   R2_gain predicted by the multi-helper expression (main-text
              equation 15) against R2_gain observed in simulation, for 300
              random architectures (100 each with 2, 5 and 10 helper traits,
              N = 50,000), with one least-squares fit across all 300
    panel B   five equivalent helper traits whose pairwise genetic and
              environmental correlation is varied from 0 to 0.95: mean of 100
              simulated datasets (N = 20,000), empirical 95% simulation interval
              (2.5th-97.5th percentiles), and the analytical prediction

The simulations were run by J.M.A.; their outputs are in reference/multi_helper/
(see the README there). This script only draws them, in the house style.
reference/multi_helper/plot_eq15_supplementary_figure_v2.R draws the same
figure in R.

Outputs (under results/figures/):
    FigS_multi_helper.{png,pdf}
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import paths  # noqa: E402
from correlated_traits.style import BLUE, MARKERS, SERIES, panel_letter, plt, save  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "reference" / "multi_helper"
HELPER_COUNTS = (2, 5, 10)


def read(name):
    with open(DATA / name, newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    paths.ensure_output_dirs()
    validation = read("eq15_validation_simulation.csv")
    summary = read("eq15_panelB_simulation_summary.csv")

    m = np.array([int(r["m"]) for r in validation])
    predicted = np.array([float(r["r2_gain_predicted"]) for r in validation])
    simulated = np.array([float(r["r2_gain_simulated"]) for r in validation])
    slope, intercept = np.polyfit(predicted, simulated, 1)
    r_squared = np.corrcoef(predicted, simulated)[0, 1] ** 2

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.8))

    # A: predicted vs simulated, one fit across all 300 architectures
    high = max(predicted.max(), simulated.max()) * 1.03
    ax_a.plot([0, high], [0, high], color="0.5", ls="--", lw=1.2, zorder=1)
    grid = np.linspace(0, high, 100)
    ax_a.plot(grid, intercept + slope * grid, color="black", lw=1.4, zorder=2)
    for count, colour, marker in zip(HELPER_COUNTS, SERIES, MARKERS):
        sel = m == count
        ax_a.plot(predicted[sel], simulated[sel], marker, ls="none", color=colour, alpha=0.7,
                  mew=0, ms=5.5, zorder=3, label="%d helpers" % count)
    ax_a.text(0.05, 0.95, "$y = %.3f + %.3fx$\n$R^2 = %.4f$" % (intercept, slope, r_squared),
              transform=ax_a.transAxes, va="top", ha="left")
    ax_a.set_xlim(0, high)
    ax_a.set_ylim(0, high)
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlabel(r"Predicted $R^2_{\mathrm{gain}}$")
    ax_a.set_ylabel(r"Simulated $R^2_{\mathrm{gain}}$")
    ax_a.legend(loc="lower right")

    # B: five equivalent helpers
    x = np.array([float(r["pairwise_helper_correlation"]) for r in summary])
    col = lambda k: np.array([float(r[k]) for r in summary])  # noqa: E731
    ax_b.fill_between(x, col("r2_gain_sim_lower95"), col("r2_gain_sim_upper95"),
                      color=BLUE, alpha=0.2, lw=0, label="95% simulation interval")
    ax_b.plot(x, col("r2_gain_analytical"), color="black", lw=1.6, label="Equation (15)")
    ax_b.plot(x, col("r2_gain_sim_mean"), "o", color=BLUE, mew=0, ms=6, label="Simulation mean")
    ax_b.set_xlim(-0.02, 0.97)
    ax_b.set_xlabel("Pairwise correlation among helper traits")
    ax_b.set_ylabel(r"$R^2_{\mathrm{gain}}$")
    ax_b.legend(loc="upper right")

    panel_letter(ax_a, "A")
    panel_letter(ax_b, "B")
    fig.tight_layout()
    out = save(fig, paths.FIGURES / "FigS_multi_helper")
    print("panel A: simulated = %.4f + %.4f x predicted, R^2 = %.4f (%d architectures)"
          % (intercept, slope, r_squared, len(m)))
    print("figure -> %s" % out)


if __name__ == "__main__":
    main()
