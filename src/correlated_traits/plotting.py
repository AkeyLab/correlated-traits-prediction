"""Shared heatmap drawing for Figures 3, 4, S1 and S2.

Colours, fonts and output resolution come from ``correlated_traits.style``.
"""

from __future__ import annotations

import numpy as np

from .style import SEQUENTIAL, panel_letter, plt, save

__all__ = ["draw_heatmap_pair", "heatmap_axes", "plt", "PANEL_SPEC"]

PANEL_SPEC = [
    ("A", r"Fix genetic correlation $\rho_g = 0.9$", r"Environmental correlation $\rho_e$"),
    ("B", r"Fix environmental correlation $\rho_e = 0.9$", r"Genetic correlation $\rho_g$"),
]


def heatmap_axes(ax, data, rho_grid, h2_h_grid, cmap, vmin, vmax, title, ylabel):
    """One heatmap panel on the (helper heritability, correlation) grid."""
    image = ax.imshow(
        data, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        extent=[h2_h_grid[0] - 0.05, h2_h_grid[-1] + 0.05, rho_grid[0], rho_grid[-1]],
    )
    ax.set_xticks(h2_h_grid)
    ax.set_xticklabels(["%.1f" % h for h in h2_h_grid])
    ax.set_yticks(np.round(np.linspace(-1, 1, 11), 1))
    ax.set_xlabel(r"Helper heritability $h_h^2$")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return image


def draw_heatmap_pair(panel_a, panel_b, rho_grid, h2_h_grid, output_path):
    """Render the two-panel heatmap used for Figures 3, 4, S1 and S2.

    Both panels share one colour scale, because the point of the pair is to
    compare them. Different *figures* must not share a scale -- the
    strong-baseline and poor-baseline scenarios differ by a factor of roughly
    fifty, and forcing a common scale renders one of them blank.
    """
    vmax = max(panel_a.max(), panel_b.max())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, data, (letter, title, ylabel) in zip(axes, (panel_a, panel_b), PANEL_SPEC):
        image = heatmap_axes(ax, data, rho_grid, h2_h_grid, SEQUENTIAL, 0, vmax, title, ylabel)
        colourbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colourbar.ax.set_title(r"$R^2_{\mathrm{gain}}$", pad=8)
        panel_letter(ax, letter)
    fig.tight_layout()
    stem = str(output_path)
    stem = stem[: -len(".png")] if stem.endswith(".png") else stem
    return save(fig, stem)
