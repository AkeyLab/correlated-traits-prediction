"""Shared plotting style for the manuscript figures.

Matplotlib is put in the non-interactive Agg backend and PDF/PS font type 42
(TrueType) so text stays selectable and editable in the vector outputs, which is
what Genome Biology's production process expects.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (import must follow backend selection)
import numpy as np  # noqa: E402

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

__all__ = ["draw_heatmap_pair", "plt"]

_PANEL_SPEC = [
    ("A", r"Fix Genetic Correlation = $0.9$", "Environmental Correlation"),
    ("B", r"Fix Environmental Correlation = $0.9$", "Genetic Correlation"),
]


def draw_heatmap_pair(panel_a, panel_b, rho_grid, h2_h_grid, output_path):
    """Render the two-panel heatmap used for Figures 3, 4, S2 and S3.

    Both panels share one colour scale, because the point of the pair is to
    compare them. Note that different *figures* must not share a scale -- the
    strong-baseline and poor-baseline scenarios differ by a factor of roughly
    fifty, and forcing a common scale renders one of them blank.
    """
    vmax = max(panel_a.max(), panel_b.max())
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    for ax, data, (letter, subtitle, ylabel) in zip(axes, (panel_a, panel_b), _PANEL_SPEC):
        image = ax.imshow(
            data, origin="lower", aspect="auto", cmap="BuPu", vmin=0, vmax=vmax,
            extent=[h2_h_grid[0] - 0.05, h2_h_grid[-1] + 0.05, rho_grid[0], rho_grid[-1]],
        )
        ax.set_xticks(h2_h_grid)
        ax.set_yticks(np.round(np.linspace(-1, 1, 11), 1))
        ax.tick_params(labelsize=13)
        ax.set_xlabel("Helper Heritability", fontsize=17)
        ax.set_ylabel(ylabel, fontsize=17)
        ax.set_title(subtitle, fontsize=16)
        ax.text(-0.14, 1.04, letter, transform=ax.transAxes, fontsize=26,
                va="top", ha="left")
        colourbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colourbar.ax.set_title(r"$R^2_{\mathrm{gain}}$", fontsize=18, pad=10)
        colourbar.ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path
