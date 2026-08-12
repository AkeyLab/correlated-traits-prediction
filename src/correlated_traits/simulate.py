"""Simulation of two correlated traits and a baseline polygenic predictor.

This implements the model written in the Methods section of the manuscript.

Genotypes
    Allele frequencies are drawn ``p_j ~ Uniform(0.05, 0.5)`` and genotypes
    ``X_ij ~ Binomial(2, p_j)``, then standardized column-wise. Variants are
    independent -- there is no linkage disequilibrium in this model.

Effect sizes
    Causal variants split into ``m_shared`` variants affecting both traits and
    ``m_helper_specific`` variants affecting only the helper. Three orthonormal
    vectors ``a``, ``b``, ``c`` are drawn on those blocks, and

        beta_t = a * sqrt(h2_t)
        beta_h = sqrt(h2_h) * (a * rho_g
                               + b * sqrt((1 - rho_g^2) * (1 - lam))
                               + c * sqrt((1 - rho_g^2) * lam))

    so the genetic correlation between the traits is ``rho_g`` by construction
    and ``lam`` is the fraction of the helper's non-shared genetic variance that
    sits on helper-specific variants.

Environment
    Environmental noise is bivariate normal with correlation ``rho_e``, scaled so
    each trait has unit total variance.

Baseline predictor
    Rather than running a genome-wide association study inside the simulation,
    the baseline is synthesized directly (Equation 13 of the manuscript):

        B = sqrt(alpha) * standardize(G_t) + sqrt(1 - alpha) * eta

    which gives ``R2_base = alpha * h2_t`` in expectation. ``alpha`` is therefore
    a dial for baseline quality: 0.1 is the poor baseline of Figures 3 and 4, and
    0.9 the strong baseline of Figures S2 and S3.

Reproducibility
    Every function takes an explicit ``numpy.random.Generator``. Results are
    bit-for-bit reproducible only if random numbers are consumed in the same
    order, so the figure drivers deliberately reuse one generator across panels
    in a fixed sequence. ``scripts/verify_reproducibility.py`` checks the
    regenerated arrays against the published ones.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .theory import joint_r2

__all__ = ["SimulationConfig", "pearson", "draw_basis", "draw_genotypes", "simulate_cell",
           "relative_gain", "heatmap_panel"]


@dataclass(frozen=True)
class SimulationConfig:
    """Fixed design shared by the heatmap figures (Figures 3, 4, S2, S3)."""

    n_individuals: int = 20_000
    m_shared: int = 100
    m_helper_specific: int = 10
    helper_specific_fraction: float = 0.20   # lam in the Methods
    n_replicates: int = 30

    @property
    def m_total(self) -> int:
        return self.m_shared + self.m_helper_specific


def pearson(x, y):
    """Pearson correlation of two 1-D arrays."""
    x = x - x.mean()
    y = y - y.mean()
    return float(x @ y / np.sqrt((x @ x) * (y @ y)))


def _rescale(v, target_variance):
    """Scale a vector to an exact sample variance. Scaling preserves correlations."""
    return v * np.sqrt(target_variance / v.var())


def draw_genotypes(rng, cfg: SimulationConfig):
    """Standardized genotype matrix, shape ``(n_individuals, m_total)``."""
    p = rng.uniform(0.05, 0.5, cfg.m_total)
    x = rng.binomial(2, p, size=(cfg.n_individuals, cfg.m_total)).astype(np.float64)
    x -= 2 * p
    x /= np.sqrt(2 * p * (1 - p))
    return x


def draw_basis(rng, cfg: SimulationConfig):
    """Orthonormal ``a``, ``b`` on the shared block; unit ``c`` on the helper block."""
    a = np.zeros(cfg.m_total)
    b = np.zeros(cfg.m_total)
    c = np.zeros(cfg.m_total)

    va = rng.standard_normal(cfg.m_shared)
    va /= np.linalg.norm(va)
    vb = rng.standard_normal(cfg.m_shared)
    vb -= (vb @ va) * va
    vb /= np.linalg.norm(vb)
    a[:cfg.m_shared] = va
    b[:cfg.m_shared] = vb

    vc = rng.standard_normal(cfg.m_helper_specific)
    vc /= np.linalg.norm(vc)
    c[cfg.m_shared:] = vc
    return a, b, c


def simulate_cell(rng, cfg, g_a, g_b, g_c, h2_t, h2_h, rho_g, rho_e, alpha):
    """Simulate one parameter cell and return ``(target, helper, baseline)``.

    All three are standardized. ``g_a``, ``g_b``, ``g_c`` are the genotype matrix
    projected onto the three basis vectors; they are reused across every cell of
    a replicate so that one genotype draw serves the whole grid.

    Random numbers are consumed here in a fixed order -- the two environmental
    draws, then the baseline noise. Callers that need extra quantities should
    derive them from the returned vectors rather than drawing again, so that the
    random stream stays comparable across scripts.
    """
    lam = cfg.helper_specific_fraction

    genetic_t = _rescale(g_a * np.sqrt(h2_t), h2_t)
    coef_b = np.sqrt(max((1 - rho_g ** 2) * (1 - lam), 0.0))
    coef_c = np.sqrt(max((1 - rho_g ** 2) * lam, 0.0))
    genetic_h = _rescale(np.sqrt(h2_h) * (g_a * rho_g + g_b * coef_b + g_c * coef_c), h2_h)

    z1 = rng.standard_normal(cfg.n_individuals)
    z2 = rng.standard_normal(cfg.n_individuals)
    env_t = np.sqrt(1 - h2_t) * z1
    env_h = np.sqrt(1 - h2_h) * (rho_e * z1 + np.sqrt(max(1 - rho_e ** 2, 0.0)) * z2)

    trait_t = genetic_t + env_t
    trait_h = genetic_h + env_h
    trait_t = (trait_t - trait_t.mean()) / trait_t.std()
    trait_h = (trait_h - trait_h.mean()) / trait_h.std()

    genetic_t_std = (genetic_t - genetic_t.mean()) / genetic_t.std()
    baseline = (np.sqrt(alpha) * genetic_t_std
                + np.sqrt(1 - alpha) * rng.standard_normal(cfg.n_individuals))

    return trait_t, trait_h, baseline


def relative_gain(rng, cfg, g_a, g_b, g_c, h2_t, h2_h, rho_g, rho_e, alpha):
    """Relative gain ``(R2_joint - R2_base) / R2_base`` for one parameter cell."""
    trait_t, trait_h, baseline = simulate_cell(
        rng, cfg, g_a, g_b, g_c, h2_t, h2_h, rho_g, rho_e, alpha
    )
    r_tb = pearson(trait_t, baseline)
    r_th = pearson(trait_t, trait_h)
    r_hb = pearson(trait_h, baseline)

    r2_base = r_tb ** 2
    return (joint_r2(r_tb, r_th, r_hb) - r2_base) / r2_base


def heatmap_panel(rng, cfg, h2_t, alpha, mode, rho_grid, h2_h_grid):
    """One heatmap panel, averaged over replicates.

    ``mode='A'`` fixes the genetic correlation at 0.9 and sweeps the
    environmental correlation down the y-axis; ``mode='B'`` fixes the
    environmental correlation at 0.9 and sweeps the genetic correlation.
    """
    if mode not in ("A", "B"):
        raise ValueError("mode must be 'A' or 'B', got %r" % (mode,))

    total = np.zeros((len(rho_grid), len(h2_h_grid)))
    for _ in range(cfg.n_replicates):
        genotypes = draw_genotypes(rng, cfg)
        a, b, c = draw_basis(rng, cfg)
        g_a, g_b, g_c = genotypes @ a, genotypes @ b, genotypes @ c
        for i, rho in enumerate(rho_grid):
            for j, h2_h in enumerate(h2_h_grid):
                rho_g, rho_e = (0.9, rho) if mode == "A" else (rho, 0.9)
                total[i, j] += relative_gain(
                    rng, cfg, g_a, g_b, g_c, h2_t, h2_h, rho_g, rho_e, alpha
                )
    return total / cfg.n_replicates
