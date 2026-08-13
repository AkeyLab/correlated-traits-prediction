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

Python 3.9 or newer. The only hard dependencies are NumPy and Matplotlib.

```bash
git clone <repository-url>
cd correlated-traits-prediction
uv sync
```

That installs the package and puts a `correlated-traits` command on the path. `uv.lock`
records an exact dependency resolution, so `uv sync` gives everyone the same versions.

<details>
<summary>Without uv</summary>

```bash
python -m venv .venv && source .venv/bin/activate
pip install .          # or: pip install -e .   for an editable checkout
```

The dependencies live in `pyproject.toml`; there is no separate requirements file.
</details>

You can also skip installation entirely — `uv run ./run.py <command>` resolves everything
from the inline metadata in `run.py` on first use.

---

## Quick start

Reproduce the cheapest figure and confirm your environment matches the published numbers
exactly. Takes about half a minute:

```bash
uv run correlated-traits quick
```

`correlated-traits --help` lists every command. Each one is also available on its own, so
`correlated-traits figure2` regenerates just that figure.

Everything is written to `results/` in the current directory (git-ignored). Redirect it with
`--results-dir PATH` — valid on either side of the command — or by setting
`CORRELATED_TRAITS_RESULTS`.

---

## Reproducing each figure

| Figure | Command | What it shows | Runtime |
|---|---|---|---|
| **2** | `correlated-traits figure2` | Simulation matches the three-term decomposition, varying one determinant per panel | seconds |
| **3, 4** | `correlated-traits figure3-4` | Predictive gain across ρ and helper heritability, **poor** baseline (α = 0.1) | 30–60 min |
| **S2, S3** | `correlated-traits figureS2-S3` | The same, **strong** baseline (α = 0.9) | 30–60 min |
| **S1** | `correlated-traits figureS1` | Closed form against simulation, with Monte Carlo error bands | seconds |

Two supporting commands quantify uncertainty, since the published heatmaps report averages
without error bars:

| Command | Purpose | Runtime |
|---|---|---|
| `correlated-traits replicate-column` | Re-runs one column with a different seed; gives per-cell standard errors and z-scores, and checks the closed-form R² identity against an explicit least-squares fit | a few min |
| `correlated-traits replicate-panel` | The same over the full grid; produces the standard errors Figure S1 plots | 10–30 min |

Dependency order — Figure S1 needs both of these first:

```
correlated-traits figure3-4        ─┐
                                    ├─→ correlated-traits figureS1
correlated-traits replicate-panel  ─┘
```

`correlated-traits figures` runs everything in the correct order and finishes by verifying
the output; `correlated-traits quick` runs only the fast commands. The grouping commands
`heatmaps` and `replications` sit in between.

`correlated-traits clean --yes` deletes the results tree.

---

## Repository layout

```
src/correlated_traits/
    theory.py          Closed-form expressions for the predictive gain
    simulate.py        Two-trait simulation model and the heatmap driver
    plotting.py        Shared figure style
    paths.py           Output locations (no absolute paths anywhere)
    cli.py             The correlated-traits command
    selftest.py        Bit-for-bit check against the original code
    figures/           One module per figure; each exposes a main()
    reference/arrays/  Published simulation arrays, for verifying reproduction
run.py              Zero-install entry point for uv run
tests/              Pytest wrapper over selftest.py
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
on the *order* in which random numbers are drawn. `src/correlated_traits/selftest.py` holds a
verbatim copy of the original implementation and asserts the package reproduces it to the
last bit.

```
$ correlated-traits test
PASS  package reproduces the original implementation bit-for-bit
PASS  two-predictor R^2 closed form matches an explicit least-squares fit
```

**2. Regenerated arrays match the published ones.** The simulation arrays behind the
published figures ship inside the package, at `src/correlated_traits/reference/arrays/`.
After regenerating, compare:

```bash
correlated-traits verify
```

Exit status is 0 when everything matches to at least floating-point tolerance, so this works
in CI. Note that a bit-for-bit match is expected on the same NumPy build; a different NumPy
or BLAS can reorder reductions and produce agreement to ~1e-16 instead, which `verify`
reports as `close` and still passes.

**3. Independent replication agrees within Monte Carlo error.** A separate seed reproduces
the closed form with a mean signed deviation of +0.1% across the column, and the closed-form
R² identity agrees with an explicit least-squares fit to ~1e-14.

**One caveat that matters.** `correlated-traits figure3-4` generates *both* figures from a
single random generator, in a fixed panel order, and `figureS2-S3` does the same for its
pair. Splitting either command, or reordering its panels, changes which random numbers each
panel receives — the results stay statistically valid but no longer match the published
arrays bit-for-bit. Each module says so at the top.

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
