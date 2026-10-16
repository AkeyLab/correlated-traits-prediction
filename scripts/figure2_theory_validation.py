#!/usr/bin/env python3
"""Figure 2 -- numerical illustration and simulation-based evaluation of theory.

Top row (A-C): simulation validates the three-term decomposition.

Equation 12 says the relative gain factorises into three independent
determinants. Each panel varies one of them and holds the other two fixed:

    Panel A   vary R2_base    (fix r_th = 0.30, redundancy = 0.30)
              gain falls as the baseline gets better
    Panel B   vary r_th       (fix R2_base = 0.20, redundancy = 0.30)
              gain is a parabola in the helper-target correlation
    Panel C   vary redundancy (fix R2_base = 0.20, r_th = 0.42)
              gain falls linearly with redundancy

These panels do not simulate genotypes. The three variables (target, helper,
baseline) are drawn directly from a multivariate normal with the required
pairwise correlations, which is enough to test the decomposition and is far
cheaper. Shaded bands are 2.5th-97.5th percentiles across 150 replicates of
2,000 individuals.

Bottom row (D-F): the closed form (File S1, Supplementary Methods 2) plotted
against the genotype-based simulation. Every simulated cell is plotted against
the value the closed form predicts for it; if the derivation is right the
points lie on the identity line. All three fix rho_g = 0.9 and alpha = 0.1.

    Panel D   Figure 3A, the single column h2_h = 0.6 (h2_t = 0.9), 21 points
    Panel E   Figure 3A, the whole panel (h2_t = 0.9), 21 x 9 = 189 cells
    Panel F   Figure 4A, the whole panel (h2_t = 0.1), 21 x 9 = 189 cells

The shaded band is the 95% Monte Carlo agreement band about the identity: where
a 30-replicate average should fall if the closed form is exact. Its half-width
is 1.96 times the standard error, and because that error grows almost linearly
with the size of the gain it is smoothed by regressing the replication standard
errors on the theoretical value.

Prerequisites for the bottom row -- run these first:
    scripts/figure3_figure4_heatmaps.py            produces the simulated values
    scripts/replicate_figure3A_figure4A_panels.py  produces the standard errors

Runtime: about 10 seconds once the prerequisites exist.
Output:  figures/Fig2_validation.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import H2_HELPER_GRID, RHO_GRID, VALIDATION_SEED, paths  # noqa: E402
from correlated_traits.plotting import plt  # noqa: E402
from correlated_traits.theory import (  # noqa: E402
    gain_closed_form, gain_from_decomposition, joint_r2,
)

N_INDIVIDUALS = 2000
N_REPLICATES = 150

LABEL_SIZE = 22
TICK_SIZE = 16
LETTER_SIZE = 34


# ----------------------------------------------------------------- top row
def solve_r_hb(r2_base, r_th, redundancy):
    """Find Cor(baseline, helper) that yields the requested redundancy.

    Redundancy is defined through the partial correlation of target and helper
    given the baseline: ``redundancy = 1 - r_th|B**2 / r_th**2``. Inverting that
    for ``r_hB`` gives a quadratic; we take the root that reproduces the correct
    sign of the partial correlation.
    """
    r_tb = np.sqrt(r2_base)
    partial = np.sign(r_th) * np.sqrt(max(r_th ** 2 * (1 - redundancy), 0.0))

    quad_a = r_tb ** 2 + partial ** 2 * (1 - r_tb ** 2)
    quad_b = -2 * r_tb * r_th
    quad_c = r_th ** 2 - partial ** 2 * (1 - r_tb ** 2)
    discriminant = max(quad_b ** 2 - 4 * quad_a * quad_c, 0.0)

    for root in ((-quad_b + np.sqrt(discriminant)) / (2 * quad_a),
                 (-quad_b - np.sqrt(discriminant)) / (2 * quad_a)):
        if abs(root) < 0.999:
            implied = (r_th - r_tb * root) / np.sqrt((1 - r_tb ** 2) * (1 - root ** 2))
            if np.sign(implied) == np.sign(partial) or partial == 0:
                return root
    return np.clip((-quad_b + np.sqrt(discriminant)) / (2 * quad_a), -0.98, 0.98)


def simulate_gains(rng, r2_base, r_th, redundancy):
    """Empirical relative gain across replicates for one parameter setting."""
    r_tb = np.sqrt(r2_base)
    r_hb = solve_r_hb(r2_base, r_th, redundancy)
    correlation = np.array([[1, r_th, r_tb], [r_th, 1, r_hb], [r_tb, r_hb, 1]])
    try:
        cholesky = np.linalg.cholesky(correlation + 1e-9 * np.eye(3))
    except np.linalg.LinAlgError:
        return np.nan

    gains = np.empty(N_REPLICATES)
    for k in range(N_REPLICATES):
        draws = (cholesky @ rng.standard_normal((3, N_INDIVIDUALS))).T
        target, helper, baseline = draws[:, 0], draws[:, 1], draws[:, 2]
        emp_r_tb = np.corrcoef(target, baseline)[0, 1]
        emp_r_th = np.corrcoef(target, helper)[0, 1]
        emp_r_hb = np.corrcoef(helper, baseline)[0, 1]
        emp_r2_base = emp_r_tb ** 2
        gains[k] = (joint_r2(emp_r_tb, emp_r_th, emp_r_hb) - emp_r2_base) / emp_r2_base
    return gains


def draw_decomposition_panel(ax, xs, gains, theoretical, colour, xlabel, letter, legend_loc):
    mean = np.array([np.mean(g) for g in gains])
    lower = np.array([np.percentile(g, 2.5) for g in gains])
    upper = np.array([np.percentile(g, 97.5) for g in gains])
    ax.fill_between(xs, lower, upper, color=colour, alpha=0.25)
    ax.plot(xs, theoretical, "k--", lw=2, label="Theory")
    ax.plot(xs, mean, color=colour, lw=2.2, label="Simulation")
    ax.set_xlabel(xlabel, fontsize=LABEL_SIZE)
    ax.set_ylabel(r"$R^2_{\mathrm{gain}}$", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=17, loc=legend_loc, framealpha=0.95)
    ax.text(-0.17, 1.06, letter, transform=ax.transAxes, fontsize=LETTER_SIZE,
            va="top", ha="left")
    return mean


#: (letter, x values, colour, x label, legend position, fixed parameters)
DECOMPOSITION_PANELS = [
    ("A", np.linspace(0.12, 0.38, 20), "#2ca02c", r"$R^2_{\mathrm{base}}$", "upper right",
     lambda x: (x, 0.30, 0.30)),
    ("B", np.linspace(-0.40, 0.40, 25), "#ff9f1c", r"$r_{th}$", "upper center",
     lambda x: (0.20, x, 0.30)),
    ("C", np.linspace(0.05, 0.65, 20), "#e63946", r"$\mathcal{R}$", "upper right",
     lambda x: (0.20, 0.42, x)),
]


# -------------------------------------------------------------- bottom row
RHO_G = 0.9
ALPHA = 0.1
COLUMN_H2_HELPER = 0.6
Z95 = 1.959964

POINT_COLOUR = "#2b6cb0"
FIT_COLOUR = "#c05621"
BAND_COLOUR = "#6b8fb8"


def simple_regression(x, y):
    """Slope, intercept and slope standard error from a simple regression."""
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    sigma_squared = (residual @ residual) / (len(x) - 2)
    inverse = np.linalg.inv(design.T @ design)
    return coefficients[1], coefficients[0], np.sqrt(sigma_squared * inverse[1, 1])


def draw_closed_form_panel(ax, predicted, simulated, standard_error, title, letter,
                           point_size):
    high = max(predicted.max(), simulated.max()) * 1.06
    grid = np.linspace(0.0, high, 400)

    se_design = np.column_stack([np.ones(len(predicted)), predicted])
    se_fit, *_ = np.linalg.lstsq(se_design, standard_error, rcond=None)
    se_curve = np.clip(se_fit[0] + se_fit[1] * grid, 0.0, None)
    ax.fill_between(grid, grid - Z95 * se_curve, grid + Z95 * se_curve,
                    color=BAND_COLOUR, alpha=0.30, lw=0, zorder=1,
                    label="95% Monte Carlo band")

    slope, intercept, slope_se = simple_regression(predicted, simulated)
    ax.plot(grid, intercept + slope * grid, color=FIT_COLOUR, lw=1.8, zorder=3,
            label="Fit: slope %.3f $\\pm$ %.3f" % (slope, slope_se))
    ax.plot([0, high], [0, high], ls=(0, (6, 5)), lw=1.6, color="0.25", zorder=3.5,
            label="Identity, $y = x$")
    ax.plot(predicted, simulated, "o", ms=point_size, mfc="white", mec=POINT_COLOUR,
            mew=1.5, zorder=4, label="Simulated")

    r_squared = np.corrcoef(predicted, simulated)[0, 1] ** 2
    ax.text(0.05, 0.95, "$R^2$ = %.4f" % r_squared, transform=ax.transAxes,
            va="top", ha="left", fontsize=16,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.75", lw=0.8))

    ax.set_xlim(0, high)
    ax.set_ylim(0, high)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Theoretical $R^2_{\mathrm{gain}}$", fontsize=LABEL_SIZE)
    ax.set_ylabel(r"Simulated $R^2_{\mathrm{gain}}$", fontsize=LABEL_SIZE)
    ax.set_title(title, fontsize=17, pad=10)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=13.5, loc="lower right", framealpha=0.95)
    ax.text(-0.17, 1.06, letter, transform=ax.transAxes, fontsize=LETTER_SIZE,
            va="top", ha="left")
    return r_squared, slope, intercept, slope_se


def require(path, produced_by):
    if not path.exists():
        raise SystemExit("Missing input: %s\nRun %s first." % (path, produced_by))
    return path


def load_design(array_name, replication_name, h2_target):
    """Simulated values, replication standard errors and closed-form predictions."""
    simulated = np.load(require(paths.ARRAYS / array_name,
                                "scripts/figure3_figure4_heatmaps.py"))
    standard_error = np.load(require(paths.TABLES / replication_name,
                                     "scripts/replicate_figure3A_figure4A_panels.py"))["se"]

    expected_shape = (len(RHO_GRID), len(H2_HELPER_GRID))
    if simulated.shape != expected_shape or standard_error.shape != expected_shape:
        raise SystemExit("Unexpected input shapes: %s and %s, expected %s"
                         % (simulated.shape, standard_error.shape, expected_shape))

    grid_rho_e, grid_h2_h = np.meshgrid(RHO_GRID, H2_HELPER_GRID, indexing="ij")
    predicted = gain_closed_form(h2_target, grid_h2_h, RHO_G, grid_rho_e,
                                 ALPHA * h2_target)
    return predicted, simulated, standard_error


def main() -> None:
    paths.ensure_output_dirs()

    # Load the bottom-row inputs first, so a missing prerequisite fails fast.
    fig3a = load_design("Fig3_scenario1_highh_A.npy", "replication_fig3A_panel.npz", 0.9)
    fig4a = load_design("Fig4_scenario2_lowh_A.npy", "replication_fig4A_panel.npz", 0.1)

    fig = plt.figure(figsize=(21, 14.2))
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.22], hspace=0.30, wspace=0.30,
                            left=0.06, right=0.985, top=0.965, bottom=0.05)

    rng = np.random.default_rng(VALIDATION_SEED)
    for column, (letter, xs, colour, xlabel, legend_loc, unpack) in enumerate(
            DECOMPOSITION_PANELS):
        gains = [simulate_gains(rng, *unpack(x)) for x in xs]
        theoretical = [gain_from_decomposition(*unpack(x)) for x in xs]
        mean = draw_decomposition_panel(fig.add_subplot(grid[0, column]), xs, gains,
                                        theoretical, colour, xlabel, letter, legend_loc)
        print("panel %s  simulated gain %.3f -> %.3f" % (letter, mean[0], mean[-1]))

    j = int(np.where(np.isclose(H2_HELPER_GRID, COLUMN_H2_HELPER))[0][0])
    closed_form_panels = [
        ("D", [a[:, j] for a in fig3a], 7.5,
         r"$h_t^2 = 0.9$,  $h_h^2 = 0.6$" "\n" r"$\rho_g = 0.9$,  $\rho_e \in [-1,\, 1]$"),
        ("E", [a.ravel() for a in fig3a], 5.0,
         r"$h_t^2 = 0.9$,  $h_h^2 \in [0.1,\, 0.9]$" "\n"
         r"$\rho_g = 0.9$,  $\rho_e \in [-1,\, 1]$"),
        ("F", [a.ravel() for a in fig4a], 5.0,
         r"$h_t^2 = 0.1$,  $h_h^2 \in [0.1,\, 0.9]$" "\n"
         r"$\rho_g = 0.9$,  $\rho_e \in [-1,\, 1]$"),
    ]
    for column, (letter, (predicted, simulated, standard_error), point_size, title) in \
            enumerate(closed_form_panels):
        r_squared, slope, intercept, slope_se = draw_closed_form_panel(
            fig.add_subplot(grid[1, column]), predicted, simulated, standard_error,
            title, letter, point_size)
        print("panel %s  R^2 = %.6f, slope = %.5f +/- %.5f, intercept = %+.5f"
              % (letter, r_squared, slope, slope_se, intercept))

    for extension in ("png", "pdf"):
        output = paths.FIGURES / ("Fig2_validation.%s" % extension)
        fig.savefig(output, dpi=200, bbox_inches="tight")
        print("figure -> %s" % output.name)
    plt.close(fig)


if __name__ == "__main__":
    main()
