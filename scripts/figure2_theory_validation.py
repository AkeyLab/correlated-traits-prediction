#!/usr/bin/env python3
"""Figure 2 -- simulation validates the three-term decomposition.

Equation 12 says the relative gain factorises into three independent
determinants. Each panel varies one of them and holds the other two fixed:

    Panel A   vary R2_base    (fix r_th = 0.30, redundancy = 0.30)
              gain falls as the baseline gets better
    Panel B   vary r_th       (fix R2_base = 0.20, redundancy = 0.30)
              gain is a parabola in the helper-target correlation
    Panel C   vary redundancy (fix R2_base = 0.20, r_th = 0.42)
              gain falls linearly with redundancy

Unlike the heatmap figures, this one does not simulate genotypes. The three
variables (target, helper, baseline) are drawn directly from a multivariate
normal with the required pairwise correlations, which is enough to test the
decomposition and is far cheaper. Runtime is a few seconds.

Shaded bands are 2.5th-97.5th percentiles across 150 replicates of 2,000
individuals.

Output:  figures/Fig2_validation.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import VALIDATION_SEED, paths  # noqa: E402
from correlated_traits.plotting import plt  # noqa: E402
from correlated_traits.theory import gain_from_decomposition, joint_r2  # noqa: E402

N_INDIVIDUALS = 2000
N_REPLICATES = 150


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


def draw_panel(ax, xs, gains, theoretical, colour, xlabel, letter, legend_loc):
    mean = np.array([np.mean(g) for g in gains])
    lower = np.array([np.percentile(g, 2.5) for g in gains])
    upper = np.array([np.percentile(g, 97.5) for g in gains])
    ax.fill_between(xs, lower, upper, color=colour, alpha=0.25)
    ax.plot(xs, theoretical, "k--", lw=2, label="Theory")
    ax.plot(xs, mean, color=colour, lw=2.2, label="Simulation")
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel(r"$R^2_{\mathrm{gain}}$", fontsize=24)
    ax.tick_params(labelsize=18)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=20, loc=legend_loc, framealpha=0.95)
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes, fontsize=36, va="top", ha="left")


#: (letter, x values, colour, x label, legend position, fixed parameters)
PANELS = [
    ("A", np.linspace(0.12, 0.38, 20), "#2ca02c", r"$R^2_{\mathrm{base}}$", "upper right",
     lambda x: (x, 0.30, 0.30)),
    ("B", np.linspace(-0.40, 0.40, 25), "#ff9f1c", r"$r_{th}$", "upper center",
     lambda x: (0.20, x, 0.30)),
    ("C", np.linspace(0.05, 0.65, 20), "#e63946", r"$\mathcal{R}$", "upper right",
     lambda x: (0.20, 0.42, x)),
]


def main() -> None:
    paths.ensure_output_dirs()
    rng = np.random.default_rng(VALIDATION_SEED)
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))

    for ax, (letter, xs, colour, xlabel, legend_loc, unpack) in zip(axes, PANELS):
        gains = [simulate_gains(rng, *unpack(x)) for x in xs]
        theoretical = [gain_from_decomposition(*unpack(x)) for x in xs]
        draw_panel(ax, xs, gains, theoretical, colour, xlabel, letter, legend_loc)
        print("panel %s  simulated gain %.3f -> %.3f"
              % (letter, np.mean(gains[0]), np.mean(gains[-1])))

    fig.tight_layout()
    output = paths.FIGURES / "Fig2_validation.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("figure -> %s" % output.name)


if __name__ == "__main__":
    main()
