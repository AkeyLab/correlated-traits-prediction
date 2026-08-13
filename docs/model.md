# The simulation model

Notation follows the Methods section of the paper.

| Symbol | Code | Meaning |
|---|---|---|
| `h²_t`, `h²_h` | `h2_t`, `h2_h` | Narrow-sense heritability of target and helper trait |
| `ρ_g` | `rho_g` | Genetic correlation between the two traits |
| `ρ_e` | `rho_e` | Environmental correlation between the two traits |
| `ρ_p` | — | Resulting phenotypic correlation |
| `α` | `alpha` | Fraction of target heritability captured by the baseline |
| `R²_base` | `r2_base` | Variance in the target explained by the baseline alone, `= α h²_t` |
| `r_th` | `r_th` | Correlation between target and helper phenotypes |
| `ℛ` | `redundancy` | Fraction of helper signal already carried by the baseline |
| `λ` | `helper_specific_fraction` | Share of the helper's non-shared genetic variance on helper-specific variants |

## Design constants

Set in `SimulationConfig`, and identical across Figures 3, 4, S2 and S3.

| Parameter | Value |
|---|---|
| Individuals per replicate | 20,000 |
| Shared causal variants | 100 |
| Helper-specific causal variants | 10 |
| `λ` | 0.20 |
| Replicates per grid cell | 30 |
| Environmental / genetic correlation axis | 21 values, −1 to 1 |
| Helper heritability axis | 9 values, 0.1 to 0.9 |

Each heatmap panel is therefore 189 cells × 30 replicates × 20,000 individuals.

## Steps in one replicate

**1. Genotypes.** `p_j ~ Uniform(0.05, 0.5)`, `X_ij ~ Binomial(2, p_j)`, then standardized
column-wise. Variants are independent — there is no linkage disequilibrium, which is the
main simplification relative to real genotype data.

**2. Effect-size basis.** Three vectors are drawn once per replicate: `a` and `b` orthonormal
on the 100 shared variants, `c` a unit vector on the 10 helper-specific ones. Genetic values
are then

```
G_t = X a √h²_t
G_h = X √h²_h ( a ρ_g  +  b √((1−ρ_g²)(1−λ))  +  c √((1−ρ_g²)λ) )
```

Because `a`, `b`, `c` are orthonormal, `Cor(G_t, G_h) = ρ_g` by construction — no rejection
sampling or numerical solve is needed. Both are rescaled to exact sample variances `h²_t` and
`h²_h`; rescaling changes no correlation.

Drawing the basis **once per replicate** and reusing it across all 189 cells is what makes
the grid affordable, and it also means neighbouring cells share genetic architecture, so the
heatmaps are smoother than 189 independent simulations would be.

**3. Environment.** `ε_t = √(1−h²_t) z₁` and
`ε_h = √(1−h²_h) (ρ_e z₁ + √(1−ρ_e²) z₂)` with `z₁, z₂` independent standard normals, giving
environmental correlation `ρ_e` and unit total trait variance. Phenotypes `Y = G + ε` are
then standardized.

**4. Baseline predictor.** Synthesized rather than estimated:

```
B = √α · standardize(G_t)  +  √(1−α) · η,      η ~ N(0,1)
```

This is the deliberate simplification of the whole design. Running a genome-wide association
study inside every cell would dominate the runtime and would confound baseline *quality* with
sample size, discovery thresholds and winner's-curse effects. Synthesizing `B` makes baseline
quality a single clean dial, `α`, with `R²_base = α h²_t` in expectation.

**5. Gain.** Fits are in-sample, with no train/test split, matching the Methods. For two
standardized predictors the joint `R²` has a closed form, so no regression is fit:

```
R²_joint = (r²_tB + r²_th − 2 r_tB r_th r_hB) / (1 − r²_hB)
R²_base  = r²_tB
R²_gain  = (R²_joint − R²_base) / R²_base
```

`scripts/replicate_figure3A_column.py` checks this closed form against an explicit
least-squares fit; they agree to about 1e-14.

## Scenarios

| Figure | Scenario | `h²_t` | `α` | Interpretation |
|---|---|---|---|---|
| 3 | 1 | 0.9 | 0.1 | Highly heritable target, poor baseline |
| 4 | 2 | 0.1 | 0.1 | Weakly heritable target, poor baseline |
| S3 | 3 | 0.1 | 0.9 | Weakly heritable target, strong baseline |
| S2 | 4 | 0.9 | 0.9 | Highly heritable target, strong baseline |

**Watch the file names — they are off by one against the manuscript.** The file stems follow
the *scenario* numbering, not the manuscript's figure numbers, and the two disagree:

| File stem | Scenario | Prints in the paper as |
|---|---|---|
| `S2_scenario3_lowh_goodbase` | 3, weakly heritable target | **Figure S3** |
| `S3_scenario4_highh_goodbase` | 4, highly heritable target | **Figure S2** |

The manuscript numbers its supplementary figures by order of appearance in Additional file 1,
where the Scenario 4 panel comes first. So reproducing "Figure S2" from this repository means
running the `S3_*` outputs, and vice versa. This has caused confusion before, including in an
earlier version of this file, which asserted the rule the wrong way round.

## Known limitations

These are properties of the model, not defects, but they bound what the simulations can show.

- **No linkage disequilibrium.** Variants are independent, so nothing here speaks to tagging,
  fine-mapping, or how polygenic scores transfer across populations.
- **One helper trait at a time.** The theory is pairwise. The empirical section handles many
  helper traits through redundancy-aware selection, but the theory is not extended to the
  multivariate case — this is stated as a limitation in the paper's Discussion.
- **Synthesized baseline.** `α` is a stand-in for everything that determines how good a
  polygenic score is. Real scores fail in structured ways this cannot capture.
- **In-sample fits.** No train/test split, so these are population-level quantities, not
  out-of-sample prediction accuracies.
- **Gaussian environment.** Environmental noise is bivariate normal with constant variance.
