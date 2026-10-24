#!/usr/bin/env python3
"""Empirical test of the helper-trait theory on UK Biobank T2D (Analysis 2).

Question: is a helper trait's incremental predictive value better explained by
its association with T2D *after conditioning on the baseline predictor* than by
its marginal association with T2D?

For every candidate helper h this computes
    1. r_th      marginal correlation between helper and T2D status
    2. r_hB      correlation between helper and the baseline linear predictor B
    3. r_th|B    partial correlation of helper and T2D given B
                 (also R = 1 - r_th|B^2 / r_th^2, the redundancy index)
    4. dAUC_h    AUC(baseline + h) - AUC(baseline) in the held-out test set,
                 both models fit on exactly the same training rows

and then compares |r_th| and |r_th|B| as predictors of dAUC_h (Spearman
correlation, with a bootstrap interval over helpers for the difference), and
lists matched pairs: helpers with nearly equal |r_th| but very different |r_th|B|.

Data splits follow the manuscript: the correlation subset (20%) gives r_th,
r_hB and r_th|B; the training subset (60%) fits B and the B + h models; the test
subset (20%) gives the AUCs. B is the baseline model's linear predictor
(log-odds) from PGS, sex, age and 10 genetic PCs, fit once on the training
subset and applied to every row.

This must run inside the UK Biobank Research Analysis Platform. It reads one
participant-level table and writes only per-helper summary statistics, which
contain no individual-level values. Do not export the input table.

Input: a CSV or Parquet file with one row per participant and columns
    t2d          0/1 T2D status
    split        "corr", "train" or "test"
    pgs, sex, age, pc1..pc10
    one column per candidate helper trait (any other column name)
Missing helper values are allowed; each helper is analysed on the rows where it
is observed, and its baseline AUC is recomputed on those same rows.

Usage:
    python helper_conditional_association.py participants.parquet --out results_dir

Caveat (see the analysis plan): the theory is for a continuous outcome and R^2;
here the outcome is binary and performance is AUC, so this tests the qualitative
prediction only, not Equation 13 quantitatively.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

BASELINE = ["pgs", "sex", "age"] + ["pc%d" % k for k in range(1, 11)]
RESERVED = {"t2d", "split", "eid"} | set(BASELINE)


def standardise(train, *others):
    mean, sd = train.mean(axis=0), train.std(axis=0)
    sd[sd == 0] = 1.0
    return [(x - mean) / sd for x in (train, *others)]


def fit_logistic(x, y):
    # L2-regularised logistic regression, as for the manuscript's baseline model
    return LogisticRegression(penalty="l2", C=1.0, max_iter=2000).fit(x, y)


def partial_corr(y, h, b):
    """Correlation of y and h after regressing each on b (with intercept)."""
    design = np.column_stack([np.ones_like(b), b])
    res_y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    res_h = h - design @ np.linalg.lstsq(design, h, rcond=None)[0]
    return np.corrcoef(res_y, res_h)[0, 1]


def per_helper(df, helpers):
    train, corr, test = (df[df.split == s] for s in ("train", "corr", "test"))
    xb_tr, xb_co, xb_te = standardise(train[BASELINE].to_numpy(float),
                                      corr[BASELINE].to_numpy(float),
                                      test[BASELINE].to_numpy(float))
    baseline = fit_logistic(xb_tr, train.t2d.to_numpy())
    b_corr = baseline.decision_function(xb_co)          # B on the correlation subset
    y_corr = corr.t2d.to_numpy(float)
    r2_base = np.corrcoef(y_corr, b_corr)[0, 1] ** 2

    rows = []
    for h in helpers:
        m_co = corr[h].notna().to_numpy()
        m_tr = train[h].notna().to_numpy()
        m_te = test[h].notna().to_numpy()
        if m_co.sum() < 100 or m_tr.sum() < 100 or m_te.sum() < 100:
            continue
        hc = corr[h].to_numpy(float)[m_co]
        r_th = np.corrcoef(y_corr[m_co], hc)[0, 1]
        r_hb = np.corrcoef(hc, b_corr[m_co])[0, 1]
        r_thb = partial_corr(y_corr[m_co], hc, b_corr[m_co])

        # baseline and baseline + h refit on the same training rows, scored on the same test rows
        x_tr = np.column_stack([xb_tr[m_tr], train[h].to_numpy(float)[m_tr]])
        x_te = np.column_stack([xb_te[m_te], test[h].to_numpy(float)[m_te]])
        x_tr, x_te = standardise(x_tr, x_te)
        y_tr, y_te = train.t2d.to_numpy()[m_tr], test.t2d.to_numpy()[m_te]
        auc_b = roc_auc_score(y_te, fit_logistic(x_tr[:, :-1], y_tr).decision_function(x_te[:, :-1]))
        auc_bh = roc_auc_score(y_te, fit_logistic(x_tr, y_tr).decision_function(x_te))

        rows.append(dict(helper=h, n_corr=int(m_co.sum()), n_test=int(m_te.sum()),
                         r_th=r_th, r_hB=r_hb, r_th_given_B=r_thb,
                         redundancy=1 - r_thb ** 2 / r_th ** 2 if r_th != 0 else np.nan,
                         auc_baseline=auc_b, auc_baseline_plus_h=auc_bh, delta_auc=auc_bh - auc_b))
    return pd.DataFrame(rows), r2_base


def compare(table, n_boot=2000, seed=1):
    """Spearman correlation of |r_th| and |r_th|B| with dAUC, and a bootstrap CI of the difference."""
    a, c, d = (table[k].abs().to_numpy() if k != "delta_auc" else table[k].to_numpy()
               for k in ("r_th", "r_th_given_B", "delta_auc"))
    rho_marg, rho_cond = spearmanr(a, d)[0], spearmanr(c, d)[0]
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(d), len(d))
        diffs.append(spearmanr(c[i], d[i])[0] - spearmanr(a[i], d[i])[0])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return dict(n_helpers=len(d), spearman_marginal=rho_marg, spearman_conditional=rho_cond,
                difference=rho_cond - rho_marg, difference_ci_low=lo, difference_ci_high=hi)


def matched_pairs(table, tol=0.01, top=10):
    """Helper pairs with |r_th| within `tol` and the largest gap in |r_th|B|."""
    rows = []
    for (_, p), (_, q) in itertools.combinations(table.iterrows(), 2):
        if abs(abs(p.r_th) - abs(q.r_th)) <= tol:
            hi, lo = (p, q) if abs(p.r_th_given_B) >= abs(q.r_th_given_B) else (q, p)
            rows.append(dict(helper_high=hi.helper, helper_low=lo.helper,
                             abs_r_th_high=abs(hi.r_th), abs_r_th_low=abs(lo.r_th),
                             abs_r_thB_high=abs(hi.r_th_given_B), abs_r_thB_low=abs(lo.r_th_given_B),
                             delta_auc_high=hi.delta_auc, delta_auc_low=lo.delta_auc,
                             gap_r_thB=abs(hi.r_th_given_B) - abs(lo.r_th_given_B)))
    return pd.DataFrame(rows).sort_values("gap_r_thB", ascending=False).head(top)


def plot(table, summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, key, label, rho in zip(
            axes, ("r_th", "r_th_given_B"),
            (r"Marginal association $|r_{th}|$", r"Conditional association $|r_{th|B}|$"),
            (summary["spearman_marginal"], summary["spearman_conditional"])):
        ax.scatter(table[key].abs(), table.delta_auc, s=22, color="#5b4a9b", alpha=0.8,
                   edgecolors="white", linewidths=0.4)
        ax.axhline(0, color="0.6", linewidth=0.8)
        ax.set_xlabel(label, fontsize=14)
        ax.set_title("Spearman %.2f" % rho, fontsize=13)
        ax.grid(color="0.92", linewidth=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$\Delta$AUC (baseline + helper $-$ baseline)", fontsize=14)
    for ax, letter in zip(axes, "AB"):
        ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=22, va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table")
    ap.add_argument("--out", default="analysis2_results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(args.table) if args.table.endswith(".parquet") else pd.read_csv(args.table)
    missing = set(BASELINE + ["t2d", "split"]) - set(df.columns)
    if missing:
        raise SystemExit("missing columns: %s" % sorted(missing))
    df = df.dropna(subset=BASELINE + ["t2d"])
    helpers = [c for c in df.columns if c not in RESERVED]

    table, r2_base = per_helper(df, helpers)
    summary = compare(table)
    summary["r2_baseline_corr_subset"] = r2_base
    pairs = matched_pairs(table)

    table.to_csv(out / "per_helper.csv", index=False)
    pd.DataFrame([summary]).to_csv(out / "summary.csv", index=False)
    pairs.to_csv(out / "matched_pairs.csv", index=False)
    plot(table, summary, out / "marginal_vs_conditional.png")
    print(pd.Series(summary).to_string())
    print("\nTop matched pairs:\n", pairs.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
