#!/usr/bin/env python3
"""Supplementary figure: predictive gain from multiple helper traits.

    panel A   the closed-form multi-helper gain (Equation 15; Supplementary
              Methods 3) against independently simulated data, for models with
              m = 2, 5 and 10 helper traits over randomly drawn genetic and
              environmental covariance structures
    panel B   the closed form for five equivalent helper traits: the joint gain
              falls as the pairwise correlation among them rises, while the gain
              from any one helper is unchanged (theory only, as in the text: "we
              then used the expression to isolate the effect of correlation")

The simulation is the genotype-based model of ``correlated_traits.simulate``
extended to m helpers. Genotypes are standardized Binomial(2, p) at 110 causal
variants; the genetic values of the target and the m helpers are G = X @ beta
with beta = Q @ chol(Sigma_G).T for a random orthonormal Q, so their genetic
covariance is Sigma_G; environmental values are multivariate normal with
covariance Sigma_E; the baseline is B = sqrt(alpha) * standardize(G_t)
+ sqrt(1 - alpha) * eta. R2_base and R2_joint are in-sample R^2 of the target on
B and on (B, helpers), as in the single-helper simulations.

Outputs (under results/):
    figures/FigS_multi_helper.{png,pdf}
    tables/multi_helper_validation.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import SimulationConfig, paths  # noqa: E402
from correlated_traits.style import BLUE, GREY, MARKERS, ORANGE, SERIES, panel_letter, plt, save  # noqa: E402
from correlated_traits.simulate import draw_genotypes  # noqa: E402

SEED = 20261024
N_INDIVIDUALS = 20_000

# panel A
HELPER_COUNTS = (2, 5, 10)
N_DESIGNS = 60          # random covariance structures per helper count
MAX_THEORY_GAIN = 6.0   # designs are redrawn above this, to keep the axis readable
N_REPLICATES_A = 10

# panel B: five equivalent helpers correlated with the target only through the
# environment (rho_g = 0), so the baseline carries none of their signal and the
# pairwise phenotypic correlation among helpers is exactly the correlation that
# enters S. These values give a single-helper gain of exactly 0.35.
M_EQUIVALENT = 5
H2_T, H2_H, ALPHA_B, RHO_G_B, RHO_E_B = 0.3, 0.3, 0.42, 0.0, 0.3
HELPER_CORRELATIONS = np.linspace(0.0, 0.95, 96)

# house style: 2, 5, 10 helpers -> blue circles, orange triangles, green squares
SERIES_STYLE = {m: (c, mk) for m, c, mk in zip(HELPER_COUNTS, SERIES, MARKERS)}


# ------------------------------------------------------------------ theory
def theory_gain(h2_t, h2_helpers, corr_g, corr_e, alpha):
    """Equations supp_multi_q_corr, supp_multi_S_corr and supp_multi_gain.

    ``corr_g`` and ``corr_e`` are (m+1) x (m+1) correlation matrices with the
    target in row/column 0.
    """
    r2_base = alpha * h2_t
    d_g, d_e = np.sqrt(h2_helpers), np.sqrt(1 - h2_helpers)
    rho_g, rho_e = corr_g[0, 1:], corr_e[0, 1:]
    q = (h2_t - r2_base) / np.sqrt(h2_t) * d_g * rho_g + np.sqrt(1 - h2_t) * d_e * rho_e
    s = (np.outer(d_g, d_g) * corr_g[1:, 1:] + np.outer(d_e, d_e) * corr_e[1:, 1:]
         - r2_base / h2_t * np.outer(d_g * rho_g, d_g * rho_g))
    return float(q @ np.linalg.solve(s, q) / r2_base)


# -------------------------------------------------------------- simulation
def r_squared(y, x):
    """In-sample R^2 of y on the columns of x (with intercept)."""
    design = np.column_stack([np.ones(len(y)), x])
    fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return 1.0 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)


def standardize(v):
    return (v - v.mean(axis=0)) / v.std(axis=0)


def simulate_gain(rng, cfg, h2_t, h2_helpers, corr_g, corr_e, alpha, helper_sets):
    """Relative gains for one replicate; ``helper_sets`` lists helper index subsets to add."""
    k = len(h2_helpers) + 1
    sd_g = np.sqrt(np.concatenate([[h2_t], h2_helpers]))
    sd_e = np.sqrt(1 - sd_g ** 2)
    genotypes = draw_genotypes(rng, cfg)
    q_basis, _ = np.linalg.qr(rng.standard_normal((cfg.m_total, k)))
    genetic = genotypes @ q_basis @ np.linalg.cholesky(np.outer(sd_g, sd_g) * corr_g).T
    environment = rng.standard_normal((cfg.n_individuals, k)) @ np.linalg.cholesky(
        np.outer(sd_e, sd_e) * corr_e).T
    traits = standardize(genetic + environment)
    baseline = (np.sqrt(alpha) * standardize(genetic[:, 0])
                + np.sqrt(1 - alpha) * rng.standard_normal(cfg.n_individuals))
    target = traits[:, 0]
    r2_base = r_squared(target, baseline)
    return [(r_squared(target, np.column_stack([baseline, traits[:, 1:][:, list(idx)]])) - r2_base)
            / r2_base for idx in helper_sets]


def random_correlation(rng, k):
    """Random correlation matrix from a two-factor model plus unique variance."""
    loadings = rng.uniform(-1, 1, size=(k, 2))
    cov = loadings @ loadings.T + np.diag(rng.uniform(0.2, 1.0, k))
    d = np.sqrt(np.diag(cov))
    return cov / np.outer(d, d)


def equicorrelated(m, rho_target, rho_helpers):
    corr = np.full((m + 1, m + 1), rho_helpers)
    corr[0, :] = corr[:, 0] = rho_target
    np.fill_diagonal(corr, 1.0)
    return corr


# ------------------------------------------------------------------ figure
def main() -> None:
    paths.ensure_output_dirs()
    cfg = SimulationConfig(n_individuals=N_INDIVIDUALS)
    rng = np.random.default_rng(SEED)

    rows = []
    for m in HELPER_COUNTS:
        for design in range(N_DESIGNS):
            # screening on the closed form only selects designs; it cannot bias
            # the comparison, because the simulation is drawn afterwards
            theory = np.inf
            while theory > MAX_THEORY_GAIN:
                h2_t, alpha = rng.uniform(0.2, 0.9, 2)
                h2_helpers = rng.uniform(0.1, 0.9, m)
                corr_g, corr_e = random_correlation(rng, m + 1), random_correlation(rng, m + 1)
                theory = theory_gain(h2_t, h2_helpers, corr_g, corr_e, alpha)
            sims = [simulate_gain(rng, cfg, h2_t, h2_helpers, corr_g, corr_e, alpha,
                                  [range(m)])[0] for _ in range(N_REPLICATES_A)]
            rows.append(dict(m=m, design=design, h2_t=h2_t, alpha=alpha,
                             theory=theory, simulated=float(np.mean(sims))))
        print("panel A: m = %d done" % m)

    five_theory, one_theory = [], []
    h2_helpers = np.full(M_EQUIVALENT, H2_H)
    for c in HELPER_CORRELATIONS:
        corr_g = equicorrelated(M_EQUIVALENT, RHO_G_B, c)
        corr_e = equicorrelated(M_EQUIVALENT, RHO_E_B, c)
        five_theory.append(theory_gain(H2_T, h2_helpers, corr_g, corr_e, ALPHA_B))
        one_theory.append(theory_gain(H2_T, h2_helpers[:1], corr_g[:2, :2], corr_e[:2, :2], ALPHA_B))
    print("panel B done")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.8))

    # A: theory vs simulation
    theory = np.array([r["theory"] for r in rows])
    simulated = np.array([r["simulated"] for r in rows])
    high = max(theory.max(), simulated.max()) * 1.05
    ax_a.plot([0, high], [0, high], "k--", lw=1.4, zorder=2)
    for m in HELPER_COUNTS:
        colour, marker = SERIES_STYLE[m]
        sel = np.array([r["m"] == m for r in rows])
        ax_a.plot(theory[sel], simulated[sel], marker, ls="none", color=colour, alpha=0.7,
                  mew=0, ms=5.5, zorder=3, label="%d helpers" % m)
    r2 = np.corrcoef(theory, simulated)[0, 1] ** 2
    ax_a.text(0.95, 0.05, "$R^2$ = %.4f" % r2, transform=ax_a.transAxes, va="bottom", ha="right")
    ax_a.set_xlim(0, high)
    ax_a.set_ylim(0, high)
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlabel(r"Theoretical $R^2_{\mathrm{gain}}$")
    ax_a.set_ylabel(r"Simulated $R^2_{\mathrm{gain}}$")
    ax_a.legend(loc="upper left")

    # B: five equivalent helpers, closed form
    x = HELPER_CORRELATIONS
    ax_b.axhline(M_EQUIVALENT * one_theory[0], color=GREY, ls=":", lw=1.6,
                 label="Five $\\times$ one helper (additive)")
    ax_b.plot(x, five_theory, color=BLUE, label="Five helpers")
    ax_b.plot(x, one_theory, color=ORANGE, ls="--", label="One helper")
    ax_b.set_xlim(0, x[-1])
    ax_b.set_ylim(0, 1.1 * M_EQUIVALENT * one_theory[0])
    ax_b.set_xlabel("Pairwise correlation among helper traits")
    ax_b.set_ylabel(r"$R^2_{\mathrm{gain}}$")
    ax_b.legend(loc="center right")

    panel_letter(ax_a, "A")
    panel_letter(ax_b, "B")
    fig.tight_layout()
    save(fig, paths.FIGURES / "FigS_multi_helper")

    with open(paths.TABLES / "multi_helper_validation.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    slope = np.polyfit(theory, simulated, 1)[0]
    print("panel A: R^2 = %.4f, slope = %.3f over %d designs" % (r2, slope, len(rows)))
    print("panel B: one helper %.3f (constant); five helpers %.3f at c = 0 -> %.3f at c = %.2f"
          % (one_theory[0], five_theory[0], five_theory[-1], x[-1]))


if __name__ == "__main__":
    main()
