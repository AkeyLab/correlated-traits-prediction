#!/usr/bin/env python3
"""Companion figure: phenotypic correlation versus predictive gain.

The phenotypic correlation rho_p depends only on h2_t, h2_h, rho_g and rho_e
(not on the baseline), so one rho_p map serves both the poor-baseline heatmaps
(alpha = 0.1, Figures 3/4) and the strong-baseline heatmaps (alpha = 0.9,
Figures S1/S2) that share a target heritability.

    panels A, B   rho_p over the same grids as the heatmap panels A and B
    panel C       R2_gain against rho_p, poor baseline   (alpha = 0.1)
    panel D       R2_gain against rho_p, strong baseline (alpha = 0.9)

Matched pairs -- grid cells with essentially the same rho_p, one correlation
mostly genetic (a panel-A cell) and the other mostly environmental (a panel-B
cell) -- are numbered in every panel.

No new simulation is run: R2_gain values are the published arrays in
reference/arrays/, and rho_p is computed from theory.phenotypic_correlation.

Usage:
    python scripts/figureS_phenotypic_correlation.py            # h2_t = 0.1 (Figures 4 and S2)
    python scripts/figureS_phenotypic_correlation.py --h2t 0.9  # h2_t = 0.9 (Figures 3 and S1)

Outputs (under results/):
    figures/FigS_rho_p_h2t_{0.1,0.9}.{png,pdf}
    tables/matched_pairs_h2t_{0.1,0.9}.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import H2_HELPER_GRID, RHO_GRID, paths, theory  # noqa: E402
from correlated_traits.plotting import PANEL_SPEC, heatmap_axes  # noqa: E402
from correlated_traits.style import BLUE, DIVERGING, MARKERS, ORANGE, panel_letter, plt, save  # noqa: E402

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "arrays"

# Array stems and baseline fractions for each target heritability.
DESIGNS = {
    0.1: {"poor": ("Fig4_scenario2_lowh", 0.1, "Fig. 4"),
          "strong": ("S2_scenario3_lowh_goodbase", 0.9, "Supplementary Fig. 2")},
    0.9: {"poor": ("Fig3_scenario1_highh", 0.1, "Fig. 3"),
          "strong": ("S3_scenario4_highh_goodbase", 0.9, "Supplementary Fig. 1")},
}

# Matched pairs as (panel, rho_g, rho_e, h2_h): first cell mostly environmental,
# second mostly genetic. Chosen from the grid with all correlations >= 0,
# |difference in rho_p| < 0.01 and the largest strong-baseline gain ratios.
PAIRS = {
    0.1: [(("B", 0.0, 0.9, 0.9), ("A", 0.9, 0.0, 0.9)),
          (("B", 0.0, 0.9, 0.8), ("A", 0.9, 0.3, 0.8)),
          (("B", 0.0, 0.9, 0.7), ("A", 0.9, 0.4, 0.6))],
    0.9: [(("B", 0.0, 0.9, 0.1), ("A", 0.9, 0.0, 0.1)),
          (("B", 0.4, 0.9, 0.1), ("A", 0.9, 0.0, 0.2)),
          (("B", 0.5, 0.9, 0.2), ("A", 0.9, 0.0, 0.3))],
}

# panel A grid -> blue circles, panel B grid -> orange triangles (house order)
SERIES = {"A": (BLUE, MARKERS[0], r"$\rho_g = 0.9$, $\rho_e$ varied (grid of A)"),
          "B": (ORANGE, MARKERS[1], r"$\rho_e = 0.9$, $\rho_g$ varied (grid of B)")}


def cell_index(panel, rho_g, rho_e):
    """Row of the heatmap array for a cell (rows sweep the non-fixed correlation)."""
    swept = rho_e if panel == "A" else rho_g
    return int(np.argmin(np.abs(RHO_GRID - swept)))


def correlations(panel, rho):
    return (0.9, rho) if panel == "A" else (rho, 0.9)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--h2t", type=float, default=0.1, choices=sorted(DESIGNS))
    h2_t = parser.parse_args().h2t
    paths.ensure_output_dirs()

    rho, h2_h = np.meshgrid(RHO_GRID, H2_HELPER_GRID, indexing="ij")
    rho_p, genetic_part, gain = {}, {}, {}
    for panel in "AB":
        rho_g, rho_e = correlations(panel, rho)
        rho_p[panel] = theory.phenotypic_correlation(h2_t, h2_h, rho_g, rho_e)
        genetic_part[panel] = rho_g * np.sqrt(h2_t * h2_h)
        for level, (stem, _alpha, _label) in DESIGNS[h2_t].items():
            gain[level, panel] = np.load(REFERENCE / ("%s_%s.npy" % (stem, panel)))

    from matplotlib import patheffects

    fig, axes = plt.subplots(2, 2, figsize=(11, 9.6))
    halo = [patheffects.withStroke(linewidth=2.5, foreground="white")]

    # A, B: rho_p on the heatmap grids, same layout as Figures 3, 4, S1, S2
    for ax, (letter, title, ylabel) in zip(axes[0], PANEL_SPEC):
        image = heatmap_axes(ax, rho_p[letter], RHO_GRID, H2_HELPER_GRID, DIVERGING, -1, 1,
                             title, ylabel)
        contours = ax.contour(H2_HELPER_GRID, RHO_GRID, rho_p[letter],
                              levels=np.arange(-0.8, 0.81, 0.2), colors="0.35",
                              linewidths=0.6, linestyles="solid")
        ax.clabel(contours, fmt=lambda v: "%.1f" % (v if abs(v) > 1e-9 else 0.0), fontsize=8)
        colourbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colourbar.ax.set_title(r"$\rho_p$", pad=8)
    for number, pair in enumerate(PAIRS[h2_t], start=1):
        for panel, rho_g, rho_e, hh in pair:
            axes[0]["AB".index(panel)].text(hh, rho_e if panel == "A" else rho_g, str(number),
                                            ha="center", va="center", fontsize=12,
                                            fontweight="bold", path_effects=halo, zorder=4)

    # C, D: R2_gain against rho_p for every grid cell
    rows = []
    for ax, level in zip(axes[1], ("poor", "strong")):
        stem, alpha, _label = DESIGNS[h2_t][level]
        top = max(gain[level, p].max() for p in "AB")
        for panel in "AB":
            colour, marker, name = SERIES[panel]
            ax.plot(rho_p[panel].ravel(), gain[level, panel].ravel(), marker, ls="none",
                    color=colour, alpha=0.6, mew=0, ms=4.5, label=name, zorder=2)
        for number, pair in enumerate(PAIRS[h2_t], start=1):
            xy = []
            for panel, rho_g, rho_e, hh in pair:
                i = cell_index(panel, rho_g, rho_e)
                j = int(np.argmin(np.abs(H2_HELPER_GRID - hh)))
                xy.append((rho_p[panel][i, j], gain[level, panel][i, j]))
                colour, marker, _ = SERIES[panel]
                ax.plot(*xy[-1], marker, color=colour, ms=8, mec="black", mew=0.8, zorder=4)
                if level == "strong":
                    rows.append(dict(pair=number, source="environmental" if panel == "B" else "genetic",
                                     panel=panel, rho_g=rho_g, rho_e=rho_e, h2_h=hh,
                                     rho_p=round(float(rho_p[panel][i, j]), 3),
                                     genetic_part=round(float(genetic_part[panel][i, j]), 3),
                                     gain_poor_baseline=round(float(gain["poor", panel][i, j]), 3),
                                     gain_strong_baseline=round(float(gain["strong", panel][i, j]), 3),
                                     gain_theory_strong=round(float(theory.gain_closed_form(
                                         h2_t, hh, rho_g, rho_e,
                                         DESIGNS[h2_t]["strong"][1] * h2_t)), 3)))
            ax.plot([xy[0][0], xy[1][0]], [xy[0][1], xy[1][1]], "k--", lw=1.0, zorder=3)
            upper = max(xy, key=lambda q: q[1])
            ax.annotate(str(number), upper, xytext=(-13, 3), textcoords="offset points",
                        fontsize=12, fontweight="bold", path_effects=halo, zorder=5)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.02 * top, 1.05 * top)
        ax.set_xlabel(r"Phenotypic correlation $\rho_p$")
        ax.set_ylabel(r"Simulated $R^2_{\mathrm{gain}}$")
        ax.set_title(r"$h_t^2 = %.1f$,  $\alpha = %.1f$" % (h2_t, alpha))
        ax.legend(loc="upper center")

    for ax, letter in zip(axes.ravel(), "ABCD"):
        panel_letter(ax, letter)
    fig.tight_layout()
    stem = "FigS_rho_p_h2t_%.1f" % h2_t
    save(fig, paths.FIGURES / stem)

    table = paths.TABLES / ("matched_pairs_h2t_%.1f.csv" % h2_t)
    with open(table, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("figure -> %s.{png,pdf}\ntable  -> %s" % (stem, table.name))
    for row in rows:
        print("  pair %(pair)d %(source)-13s panel %(panel)s rho_g=%(rho_g)+.1f rho_e=%(rho_e)+.1f "
              "h2_h=%(h2_h).1f rho_p=%(rho_p).3f  gain poor=%(gain_poor_baseline)7.3f "
              "strong=%(gain_strong_baseline)6.3f (theory %(gain_theory_strong)6.3f)" % row)


if __name__ == "__main__":
    main()
