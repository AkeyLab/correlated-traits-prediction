# Correlated traits improve risk prediction

Theory and simulation code for **"A Principled Framework for Using Correlated Traits to
Improve Risk Prediction"** (Zhang, Bierman & Akey).

Polygenic scores usually predict one trait from genetics alone. This paper asks when adding
a second, correlated **helper trait** to the model actually helps — and shows that the
answer is not determined by how strongly the two traits correlate, but by *why* they
correlate. This repository contains the derivations in executable form, the simulations that
validate them, and the code that produces every theoretical figure in the paper.

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

We've created an [interactive webpage](https://akeylab.github.io/correlated-traits-prediction)
to help you build intuition.

---

## Scope — read this before you look for something that is not here

This repository covers the **theoretical and simulation** results: Figures 2, 3, 4, S1, S2
and S3, plus the derivation checks.

It does **not** contain the UK Biobank type 2 diabetes analysis (Figure 5, and the
AUC-ROC values 0.677 / 0.907 / 0.889 quoted in the abstract). That analysis runs on the UK
Biobank Research Analysis Platform against controlled-access individual-level data, and lives
inside that environment. See [`docs/uk_biobank_analysis.md`](docs/uk_biobank_analysis.md) for
exactly what is missing, how to request the underlying data, and what releasing that code
would involve.

---

## Installation

This analysis ships as a Python package, `correlated_traits`, which installs one
command-line program: **`correlated-traits`**. Each subcommand regenerates one piece of the
paper such as `correlated-traits figure2` writes Figure 2.

> **Note — pre-merge URLs.** The commands below install from the `repackage` branch, because
> that is where the packaged CLI currently lives. Once it is merged, drop the `@repackage`
> suffix from every URL on this page so they track the default branch.

### With uv

If you use [uv](https://docs.astral.sh/uv/) then you can install with:

```bash
uv tool install git+https://github.com/AkeyLab/correlated-traits-prediction@repackage
correlated-traits --help
```

<details>
<summary>Without uv</summary>

Install into a virtual environment:

```bash
python -m venv correlated-traits-venv
source correlated-traits-venv/bin/activate
pip install git+https://github.com/AkeyLab/correlated-traits-prediction@repackage
correlated-traits --help
```

</details>

---

## Quick start

After you have `correlated-traits --help` working, try the cheapest figure. It also runs the
environment check, so it confirms your NumPy produces the expected numbers before you commit
to anything longer. Takes about five seconds:

```bash
correlated-traits quick
```

Everything is written to `results/` in the current directory.

---

## Reproducing each figure

| Figure | Command | What it shows | Runtime |
|---|---|---|---|
| **2** | `correlated-traits figure2` | Simulation matches the three-term decomposition, varying one determinant per panel | 4 s |
| **3, 4** | `correlated-traits figure3-4` | Predictive gain across ρ and helper heritability, **poor** baseline (α = 0.1) | 50 s |
| **S2, S3** | `correlated-traits figureS2-S3` | The same, **strong** baseline (α = 0.9) | 50 s |
| **S1** | `correlated-traits figureS1` | Closed form against simulation, with Monte Carlo error bands | 2 s |

Two supporting commands quantify uncertainty, since the published heatmaps report averages
without error bars:

| Command | Purpose | Runtime |
|---|---|---|
| `correlated-traits replicate-column` | Re-runs one column with a different seed; gives per-cell standard errors and z-scores, and checks the closed-form R² identity against an explicit least-squares fit | 6 s |
| `correlated-traits replicate-panel` | The same over the full grid; produces the standard errors Figure S1 plots | 13 s |

Dependency order — Figure S1 needs both of these first:

```
correlated-traits figure3-4        ─┐
                                    ├─→ correlated-traits figureS1
correlated-traits replicate-panel  ─┘
```

`correlated-traits figures` runs everything in the correct order and finishes by verifying
the output, and takes about two minutes end to end. `correlated-traits quick` runs only the
fastest pair; the grouping commands `heatmaps` and `replications` sit in between.

Runtimes above were measured on one core of a 2026 cluster node. They are far shorter than
the estimates this repository previously carried, which were never measured.

`correlated-traits clean --yes` deletes the results tree.

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

Results are bit-for-bit reproducible, and the repository checks this rather than asserting it.

**1. Regenerated arrays match the published ones.** The simulation arrays behind the four
heatmap figures ship inside the package. After regenerating, compare against them:

```bash
correlated-traits verify
```

Every array is reported as `identical`, `close`, or `DIFFERS`. Exit status is 0 when all of
them match to at least floating-point tolerance and 1 otherwise, so this works in CI.
`correlated-traits figures` ends by running it.

**2. A two-second check that your install is sound.** Worth running first:

```
$ correlated-traits test
PASS  simulation reproduces the reference implementation bit-for-bit
PASS  two-predictor R^2 closed form matches an explicit least-squares fit
PASS  8 published reference arrays ship with the package
```

The first line runs the simulation at reduced scale against a frozen reference implementation
embedded in the package, confirming your NumPy produces the expected random stream. The
second independently confirms that the two-predictor R² closed form the theory relies on
agrees with an explicit least-squares fit. The third confirms the published arrays that
`verify` compares against are actually present, so an incomplete install surfaces
immediately rather than at the end of a full run.

**3. Independent replication agrees within Monte Carlo error.** A separate seed reproduces
the closed form with a mean signed deviation of +0.1% across the column, and the closed-form
R² identity agrees with an explicit least-squares fit to ~1e-14.

**Which versions give an exact match.** The published figures were produced with
**Python 3.10.9, NumPy 1.23.1 and Matplotlib 3.7.2** on Linux. Use those if you want `verify`
to report `identical` for all eight arrays.

Newer NumPy should reproduce them too, and usually reports `identical` as well. `default_rng`
guarantees a stable stream for the three distributions used here — `uniform`, `binomial` and
`standard_normal` — so the draws themselves do not move between versions. What can move is
the last bit of a sum, when a different NumPy or BLAS build reorders a reduction; `verify`
then reports `close` instead and still exits 0. On NumPy 2.x the observed disagreement is at
most 5e-14 in absolute terms, on gain values ranging up to about 100 — relative agreement
near machine precision, which affects no conclusion in the paper. `verify` fails on anything
larger than 1e-10.

NumPy 1.17 is a hard floor rather than a formality: below it `default_rng` does not exist, and
the legacy `RandomState` API draws entirely different numbers.

**One caveat that matters.** `correlated-traits figure3-4` generates *both* figures from a
single random generator, in a fixed panel order, and `figureS2-S3` does the same for its
pair. Splitting either command, or reordering its panels, changes which random numbers each
panel receives — the results stay statistically valid but no longer match the published
arrays. Each module says so at the top.

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
