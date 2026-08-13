# Provenance

This repository was assembled from a working research tree that had accumulated since 2023.
This file records exactly what was brought across, what was changed, and what was left
behind, so that nothing looks like it appeared from nowhere.

## Source of each file

| This repository | Came from |
|---|---|
| `src/correlated_traits/simulate.py` | The simulation body shared by `regenerate_fig34_methods.py` and `regenerate_S2S3_methods.py` |
| `src/correlated_traits/theory.py` | The closed forms duplicated across four scripts |
| `src/correlated_traits/plotting.py` | `replot_fig34.py`, which produced the published heatmap styling |
| `src/correlated_traits/figures/figure2_theory_validation.py` | `regenerate_fig2_validation.py` |
| `src/correlated_traits/figures/figure3_figure4_heatmaps.py` | `regenerate_fig34_methods.py` |
| `src/correlated_traits/figures/figureS2_figureS3_heatmaps.py` | `regenerate_S2S3_methods.py` |
| `src/correlated_traits/figures/figureS1_theory_vs_simulation.py` | `make_fig_theory_vs_sim.py` |
| `src/correlated_traits/figures/replicate_figure3A_column.py` | `replicate_fig3A_column.py` |
| `src/correlated_traits/figures/replicate_figure3A_panel.py` | `replicate_fig3A_panel.py` |
| `src/correlated_traits/reference/arrays/*.npy` | `fig34_methods_out/*.npy`, the arrays the published figures were drawn from |

## What changed

Nothing that affects a number. The changes are structural:

1. **Deduplication.** The closed-form gain expression appeared verbatim in four scripts and
   the simulation body in three. Each now has one definition.
2. **Absolute paths removed.** Scripts previously wrote to a hard-coded home directory. All
   output now goes to `results/`, overridable with `CORRELATED_TRAITS_RESULTS`.
3. **Module-level side effects removed.** Several scripts ran on import, which made them
   impossible to reuse. Each now exposes a `main()` that the `correlated-traits` CLI calls.
4. **Global random generators replaced.** Generators are created in `main()` and passed
   explicitly, so the random stream is visible rather than implicit.
5. **Names spelled out.** `gain_cell` → `relative_gain`, `ols_r2` →
   `ordinary_least_squares_r2`, `H_T2` → `h2_t`, and so on.
6. **Prerequisites checked.** Figure S1 previously failed with a bare `FileNotFoundError` if
   run out of order; it now names the command you need to run first.

### Verification that the numbers did not move

The refactor is checked three ways, all reproducible:

- `src/correlated_traits/selftest.py` holds a verbatim copy of the original implementation and
  asserts the package reproduces it **bit-for-bit** on the same seed.
- Regenerating Figures 3 and 4 reproduced the published `.npy` arrays; check with
  `correlated-traits verify`.
- Re-running the column replication reproduced its original CSV byte-for-byte.
- Figure 2's reported gain values match the original run log exactly
  (0.479 → 0.104, 0.458, 0.669 → 0.255).

## What was deliberately left out

The source tree contains roughly 88 notebooks and several hundred scripts spanning three
years of exploration. Almost all of it is superseded. Only code that produces a figure or a
check appearing in the submitted manuscript was brought across.

Specifically not included, and why:

| Left out | Reason |
|---|---|
| `20230*`–`20260*` exploratory notebooks | Superseded working notebooks. Earlier simulation schemes, abandoned parameterizations, and figure drafts that did not reach the paper. |
| `simulation_class.py`, `simulation_class_sort_h2.py` and the `sim2_*`, `sim_part*` families | An earlier simulation framework. The published Figures 3, 4, S2 and S3 were regenerated from the Methods-faithful code, deliberately replacing this variant-selection approach. |
| `run_twoModels*.py`, `sim_all_parameters.py` | Meta-modelling that predicts the simulated gain from its own input parameters. This produced an XGBoost figure that is **not** in the current manuscript. |
| `optuna_gb.py`, `example_optuna.py` | Hyperparameter search experiments, not used in the paper. |
| `make_combined_pdf.py`, `plot_appendix2_scatter.py` | One-off assembly scripts for a specific internal review document. `plot_appendix2_scatter.py` largely duplicated the Figure S1 script. |
| `verify_appendix2_fig3A.py` | An earlier version of the Figure S1 check, superseded by the two replication scripts, which quantify uncertainty rather than only comparing point values. |
| `correlated_traits_venv/` | A checked-in virtual environment. Replaced by the dependencies declared in `pyproject.toml`. |
| Everything under the UK Biobank analysis | Controlled-access data. See [`uk_biobank_analysis.md`](uk_biobank_analysis.md). |

None of the excluded files were deleted from the original tree; this repository is a
curated copy, not a move.

## Environment the published figures were produced in

Python 3.10.9, NumPy 1.23.1, Matplotlib 3.7.2, on Linux.

NumPy's `default_rng` stream is stable across versions for the distributions used here
(`uniform`, `binomial`, `standard_normal`), so a bit-for-bit match is expected on other
recent NumPy versions too. If `correlated-traits verify` reports differences at the 1e-15
level, that is a BLAS reduction-order difference and does not affect any conclusion.
