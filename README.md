# Correlated traits improve risk prediction

Theory and simulation code for **"A Principled Framework for Using Correlated Traits to
Improve Risk Prediction"** (Zhang, Bierman & Akey).

Polygenic scores usually predict one trait from genetics alone. This paper asks when adding
a second, correlated **helper trait** to the model actually helps — and shows that the
answer is not determined by how strongly the two traits correlate, but by *why* they
correlate. This repository contains the derivations in executable form, the simulations that
validate them, and the scripts that produce every theoretical figure in the paper.

---

## The result in one equation

The relative gain in variance explained from adding a helper trait to a baseline predictor
factors into three independent pieces:

```
R²_gain  =  (1 − R²_base) / R²_base  ×  r²_th  ×  (1 − ℛ)
             ─────────────────────      ─────      ───────
             room left to improve       helper     signal not already
                                        signal     in the baseline
```

Written instead in terms of biology — heritabilities and the genetic and environmental
correlations — the same quantity becomes

```
              [ ρ_g √(h²_h/h²_t) (h²_t − R²_base)  +  ρ_e √((1−h²_t)(1−h²_h)) ]²
R²_gain  =   ──────────────────────────────────────────────────────────────────
                        R²_base [ 1 − ρ_g² R²_base h²_h / h²_t ]
```

Both forms live in [`src/correlated_traits/theory.py`](src/correlated_traits/theory.py).
The second is derived in Supplementary Methods 2 of the paper and validated against
simulation in Figure S1, where it recovers the simulated values with R² > 0.999.

The counter-intuitive consequence, and the reason the paper exists: two helper traits with
*identical* phenotypic correlation to the target can deliver completely different gains. When
the target trait is weakly heritable, a helper trait with **low** heritability is often the
better choice, because its value lies in environmental signal the polygenic baseline cannot
already contain.

---

## Scope — read this before you look for something that is not here

This repository covers the **theoretical and simulation** results: Figures 2, 3, 4, S1, S2
and S3, plus the derivation checks.

It does **not** contain the UK Biobank type 2 diabetes analysis (Figure 5, and the
AUC-ROC values 0.677 / 0.907 / 0.889 quoted in the abstract). That analysis runs against
controlled-access individual-level data and is not part of this codebase. See
[`docs/uk_biobank_analysis.md`](docs/uk_biobank_analysis.md) for what would be needed to add
it and how to request the underlying data.

---

## Installation

Python 3.9 or newer. The only hard dependencies are NumPy and Matplotlib.

```bash
git clone <repository-url>
cd correlated-traits-prediction
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

No installation step is required — the scripts add `src/` to the path themselves. To import
the package from elsewhere, `pip install -e .`

---

## Quick start

Reproduce the cheapest figure, which takes a few seconds and exercises the whole theory path:

```bash
python scripts/figure2_theory_validation.py
```

Then confirm your environment reproduces the published numbers exactly:

```bash
python tests/test_equivalence.py
```

Everything is written to `results/` (git-ignored). Set `CORRELATED_TRAITS_RESULTS` to
redirect output elsewhere.

---

## Reproducing each figure

| Figure | Script | What it shows | Runtime |
|---|---|---|---|
| **2** | `scripts/figure2_theory_validation.py` | Simulation matches the three-term decomposition, varying one determinant per panel | ~10 s |
| **3, 4** | `scripts/figure3_figure4_heatmaps.py` | Predictive gain across ρ and helper heritability, **poor** baseline (α = 0.1) | ~25 min |
| **S2, S3** | `scripts/figureS2_figureS3_heatmaps.py` | The same, **strong** baseline (α = 0.9) | ~25 min |
| **S1** | `scripts/figureS1_theory_vs_simulation.py` | Closed form against simulation, with Monte Carlo error bands | ~5 s |

Two supporting scripts quantify uncertainty, since the published heatmaps report averages
without error bars:

| Script | Purpose | Runtime |
|---|---|---|
| `scripts/replicate_figure3A_column.py` | Re-runs one column with a different seed; gives per-cell standard errors and z-scores, and checks the closed-form R² identity against an explicit least-squares fit | ~2 min |
| `scripts/replicate_figure3A_panel.py` | The same over the full grid; produces the standard errors Figure S1 plots | ~20 min |

Dependency order — Figure S1 needs both of these first:

```
figure3_figure4_heatmaps.py ─┐
                             ├─→ figureS1_theory_vs_simulation.py
replicate_figure3A_panel.py ─┘
```

`make figures` runs everything in the correct order; `make quick` runs only the fast targets.

---

## Repository layout

```
src/correlated_traits/
    theory.py       Closed-form expressions for the predictive gain
    simulate.py     Two-trait simulation model and the heatmap driver
    plotting.py     Shared figure style
    paths.py        Output locations (no absolute paths anywhere)
scripts/            One script per figure; each is runnable on its own
tests/              Bit-for-bit regression test against the original code
reference/arrays/   Published simulation arrays, for verifying reproduction
docs/               Model notes, figure map, UK Biobank analysis status
results/            All output lands here (git-ignored)
```

---

## The simulation model

Each replicate simulates 20,000 individuals.

**Genotypes.** Allele frequencies `p_j ~ Uniform(0.05, 0.5)`, genotypes
`X_ij ~ Binomial(2, p_j)`, standardized. Variants are independent — there is no linkage
disequilibrium.

**Effect sizes.** 100 causal variants shared between the traits and 10 helper-specific ones.
Three orthonormal vectors `a`, `b`, `c` are drawn across those blocks so that the genetic
correlation equals `ρ_g` by construction:

```
β_t = a √h²_t
β_h = √h²_h ( a ρ_g  +  b √((1−ρ_g²)(1−λ))  +  c √((1−ρ_g²)λ) )
```

**Environment.** Bivariate normal with correlation `ρ_e`, scaled so each trait has unit
variance.

**Baseline predictor.** Rather than running a genome-wide association study inside the
simulation, the baseline is synthesized directly:

```
B = √α · standardize(G_t)  +  √(1−α) · η
```

so `R²_base = α h²_t` in expectation. This makes α a clean dial for baseline quality: **0.1**
is the poor baseline of Figures 3 and 4, **0.9** the strong baseline of Figures S2 and S3.

Full detail in [`docs/model.md`](docs/model.md).

---

## Reproducibility

Results are bit-for-bit reproducible, and the repository checks this three ways rather than
asserting it.

**1. The package matches the original code exactly.** The published figures came from
standalone scripts that interleaved simulation, parameters and plotting in single files. This
repository factors that apart, which is risky: results depend not only on the arithmetic but
on the *order* in which random numbers are drawn. `tests/test_equivalence.py` holds a verbatim
copy of the original implementation and asserts the package reproduces it to the last bit.

```
$ python tests/test_equivalence.py
PASS  package reproduces the original implementation bit-for-bit
PASS  two-predictor R^2 closed form matches an explicit least-squares fit
```

**2. Regenerated arrays match the published ones.** `reference/arrays/` holds the simulation
arrays behind the published figures. After regenerating, compare:

```bash
python scripts/verify_reproducibility.py
```

**3. Independent replication agrees within Monte Carlo error.** A separate seed reproduces
the closed form with a mean signed deviation of +0.1% across the column, and the closed-form
R² identity agrees with an explicit least-squares fit to ~1e-14.

**One caveat that matters.** `figure3_figure4_heatmaps.py` generates *both* figures from a
single random generator, in a fixed panel order, and `figureS2_figureS3_heatmaps.py` does the
same for its pair. Splitting either script, or reordering its panels, changes which random
numbers each panel receives — the results stay statistically valid but no longer match the
published arrays bit-for-bit. Each script says so at the top.

---

## Citation

```bibtex
@article{zhang2026correlated,
  title   = {A Principled Framework for Using Correlated Traits to Improve Risk Prediction},
  author  = {Zhang, Kaiqian and Bierman, Robert and Akey, Joshua M.},
  journal = {Genome Biology},
  year    = {2026},
  note    = {Under review}
}
```

See [`CITATION.cff`](CITATION.cff). Update both when the paper is accepted.

---

## License

MIT — see [`LICENSE`](LICENSE).
