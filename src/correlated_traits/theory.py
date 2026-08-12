"""Closed-form expressions for the relative gain in predictive accuracy.

Two equivalent forms of the same quantity appear in the manuscript.

**Decomposition form** (Equation 12, first line). The gain is written in terms of
three *observable* quantities -- how good the baseline already is, how much signal
the helper trait carries, and how much of that signal the baseline already has:

    R2_gain = (1 - R2_base) / R2_base * r_th**2 * (1 - redundancy)

**Closed form** (Equation 12, second line; derived in Supplementary Methods 2).
The same gain written in terms of *biological* quantities -- heritabilities and
the genetic and environmental correlations:

                [rho_g * sqrt(h2_h / h2_t) * (h2_t - R2_base)
                     + rho_e * sqrt((1 - h2_t) * (1 - h2_h))]**2
    R2_gain = ---------------------------------------------------
                R2_base * [1 - rho_g**2 * R2_base * h2_h / h2_t]

Symbols
-------
h2_t, h2_h   Narrow-sense heritability of the target and helper trait.
rho_g        Genetic correlation between target and helper.
rho_e        Environmental correlation between target and helper.
R2_base      Variance in the target explained by the baseline predictor alone.
alpha        Fraction of the target heritability the baseline captures, so that
             R2_base = alpha * h2_t.
r_th         Pearson correlation between target and helper phenotypes.
redundancy   Fraction of the helper's signal already carried by the baseline.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "phenotypic_correlation",
    "gain_from_decomposition",
    "gain_closed_form",
    "joint_r2",
]


def phenotypic_correlation(h2_t, h2_h, rho_g, rho_e):
    """Population phenotypic correlation between two standardized traits.

    Genetic and environmental correlations combine weighted by the square roots
    of the heritable and non-heritable variance fractions. Note that opposite
    signs partially cancel, which is what produces the asymmetry visible in the
    heatmap figures.
    """
    return (rho_g * np.sqrt(h2_t * h2_h)
            + rho_e * np.sqrt((1.0 - h2_t) * (1.0 - h2_h)))


def gain_from_decomposition(r2_base, r_th, redundancy):
    """Relative gain from the three-term decomposition (Equation 12, first line)."""
    return (1.0 - r2_base) / r2_base * r_th ** 2 * (1.0 - redundancy)


def gain_closed_form(h2_t, h2_h, rho_g, rho_e, r2_base):
    """Relative gain in terms of heritabilities and correlations.

    This is the expression validated against simulation in Figure S1.
    """
    numerator = (rho_g * np.sqrt(h2_h / h2_t) * (h2_t - r2_base)
                 + rho_e * np.sqrt((1.0 - h2_t) * (1.0 - h2_h))) ** 2
    denominator = r2_base * (1.0 - rho_g ** 2 * r2_base * h2_h / h2_t)
    return numerator / denominator


def joint_r2(r_tb, r_th, r_hb):
    """In-sample R^2 of the target on two standardized predictors.

    For two predictors this squared multiple correlation has a closed form, so
    no regression needs to be fit. ``scripts/replicate_figure3A_column.py``
    checks this against an explicit least-squares fit and they agree to ~1e-14.
    """
    return (r_tb ** 2 + r_th ** 2 - 2.0 * r_tb * r_th * r_hb) / (1.0 - r_hb ** 2)
