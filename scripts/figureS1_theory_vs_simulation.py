#!/usr/bin/env python3
"""Figure S1 -- the closed form plotted against simulation.

Every simulated cell of Figure 3 panel A is plotted against the value the
closed form (Supplementary Methods 2) predicts for it. If the derivation is
right the points lie on the identity line.

    Panel A   the single column h2_h = 0.6, 21 points as rho_e sweeps -1 to 1
    Panel B   the whole panel, all 21 x 9 = 189 cells

The shaded band is the 95% Monte Carlo agreement band about the identity: where
a 30-replicate average should fall if the closed form is exact. Its half-width
is 1.96 times the standard error, and because that error grows almost linearly
with the size of the gain it is smoothed by regressing the replication standard
errors on the theoretical value.

Prerequisites -- run these first:
    scripts/figure3_figure4_heatmaps.py     produces the simulated values
    scripts/replicate_figure3A_panel.py     produces the standard errors

Runtime: seconds.
Output:  figures/FigureS1_theory_vs_simulation.{pdf,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import H2_HELPER_GRID, RHO_GRID, paths  # noqa: E402
from correlated_traits.plotting import plt  # noqa: E402
from correlated_traits.theory import gain_closed_form  # noqa: E402

H2_TARGET = 0.9
RHO_G = 0.9
ALPHA = 0.1
R2_BASE = ALPHA * H2_TARGET
COLUMN_H2_HELPER = 0.6
Z95 = 1.959964

POINT_COLOUR = "#2b6cb0"
FIT_COLOUR = "#c05621"
BAND_COLOUR = "#6b8fb8"

SIMULATED_ARRAY = "Fig3_scenario1_highh_A.npy"
REPLICATION_NPZ = "replication_fig3A_panel.npz"


def simple_regression(x, y):
    """Slope, intercept and slope standard error from a simple regression."""
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    sigma_squared = (residual @ residual) / (len(x) - 2)
    inverse = np.linalg.inv(design.T @ design)
    return coefficients[1], coefficients[0], np.sqrt(sigma_squared * inverse[1, 1])


def draw_panel(ax, predicted, simulated, standard_error, title, letter, point_size):
    low = 0.0
    high = max(predicted.max(), simulated.max()) * 1.06
    grid = np.linspace(low, high, 400)

    se_design = np.column_stack([np.ones(len(predicted)), predicted])
    se_fit, *_ = np.linalg.lstsq(se_design, standard_error, rcond=None)
    se_curve = np.clip(se_fit[0] + se_fit[1] * grid, 0.0, None)
    ax.fill_between(grid, grid - Z95 * se_curve, grid + Z95 * se_curve,
                    color=BAND_COLOUR, alpha=0.30, lw=0, zorder=1,
                    label="95% Monte Carlo band about $y = x$")

    slope, intercept, slope_se = simple_regression(predicted, simulated)
    ax.plot(grid, intercept + slope * grid, color=FIT_COLOUR, lw=1.6, zorder=3,
            label="Fit: slope %.4f $\\pm$ %.4f, intercept %+.4f"
                  % (slope, slope_se, intercept))
    ax.plot([low, high], [low, high], ls=(0, (6, 5)), lw=1.4, color="0.25",
            zorder=3.5, label="Identity, $y = x$")
    ax.plot(predicted, simulated, "o", ms=point_size, mfc="white",
            mec=POINT_COLOUR, mew=1.4, zorder=4, label="Simulated")

    r_squared = np.corrcoef(predicted, simulated)[0, 1] ** 2
    ax.text(0.04, 0.955, "$R^2$ = %.5f" % r_squared, transform=ax.transAxes,
            va="top", ha="left", fontsize=12,
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="0.75", lw=0.8))

    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Theoretical $R^2_{\mathrm{gain}}$ (Supplementary Methods 2)", fontsize=14)
    ax.set_ylabel(r"Simulated $R^2_{\mathrm{gain}}$ (Figure 3A)", fontsize=14)
    ax.set_title(title, fontsize=13.5, pad=10)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=9.5, loc="lower right", framealpha=0.95)
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes, fontsize=26, va="top", ha="left")
    return r_squared, slope, intercept, slope_se


def require(path, produced_by):
    if not path.exists():
        raise SystemExit(
            "Missing input: %s\nRun %s first." % (path, produced_by)
        )
    return path


def main() -> None:
    paths.ensure_output_dirs()

    simulated_panel = np.load(require(
        paths.ARRAYS / SIMULATED_ARRAY, "scripts/figure3_figure4_heatmaps.py"))
    standard_error_panel = np.load(require(
        paths.TABLES / REPLICATION_NPZ, "scripts/replicate_figure3A_panel.py"))["se"]

    expected_shape = (len(RHO_GRID), len(H2_HELPER_GRID))
    if simulated_panel.shape != expected_shape or standard_error_panel.shape != expected_shape:
        raise SystemExit("Unexpected input shapes: %s and %s, expected %s"
                         % (simulated_panel.shape, standard_error_panel.shape, expected_shape))

    grid_rho_e, grid_h2_h = np.meshgrid(RHO_GRID, H2_HELPER_GRID, indexing="ij")
    predicted_panel = gain_closed_form(H2_TARGET, grid_h2_h, RHO_G, grid_rho_e, R2_BASE)
    column = int(np.where(np.isclose(H2_HELPER_GRID, COLUMN_H2_HELPER))[0][0])

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 7.0))
    statistics = {
        "A": draw_panel(
            axes[0], predicted_panel[:, column], simulated_panel[:, column],
            standard_error_panel[:, column],
            r"$h_h^2 = 0.6$,  $h_t^2 = 0.9$,  $\rho_g = 0.9$" "\n"
            r"$\rho_e$ varied over $[-1,\, 1]$", "A", 7.0),
        "B": draw_panel(
            axes[1], predicted_panel.ravel(), simulated_panel.ravel(),
            standard_error_panel.ravel(),
            r"$h_h^2$ varied over $[0.1,\, 0.9]$,  $h_t^2 = 0.9$,  $\rho_g = 0.9$" "\n"
            r"$\rho_e$ varied over $[-1,\, 1]$", "B", 4.5),
    }

    fig.tight_layout()
    for extension in ("pdf", "png"):
        output = paths.FIGURES / ("FigureS1_theory_vs_simulation.%s" % extension)
        fig.savefig(output, dpi=200, bbox_inches="tight")
        print("wrote %s" % output)
    plt.close(fig)

    for tag, (r_squared, slope, intercept, slope_se) in statistics.items():
        print("  panel %s: R^2 = %.6f, slope = %.5f +/- %.5f, intercept = %+.5f"
              % (tag, r_squared, slope, slope_se, intercept))


if __name__ == "__main__":
    main()
