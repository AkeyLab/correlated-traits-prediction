"""House style shared by every code-generated figure (Figures 2-4, S1-S4).

The palette is matplotlib's default categorical cycle (tab10), always used in
the same order -- blue, orange, green -- with a matching marker for each, so a
series looks the same wherever it appears:

    BLUE   #1f77b4   circles
    ORANGE #ff7f0e   triangles
    GREEN  #2ca02c   squares

Magnitudes (R2_gain heatmaps) use one hue, white to BLUE to a darker BLUE.
Signed quantities (rho_p) use a diverging scale built from the same two
palette colours, ORANGE - white - BLUE, with white at zero.

Every figure is written twice: a 600 dpi PNG for the manuscript and a vector
PDF for production.
"""

from __future__ import annotations

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

__all__ = ["BLUE", "ORANGE", "GREEN", "GREY", "SERIES", "MARKERS", "SEQUENTIAL", "DIVERGING",
           "LETTER_SIZE", "apply", "panel_letter", "save", "plt"]

BLUE, ORANGE, GREEN, GREY = "#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f"
SERIES = (BLUE, ORANGE, GREEN)
MARKERS = ("o", "^", "s")

DARK_BLUE = "#0f3c5a"   # BLUE at half brightness, for the top of the sequential scale
SEQUENTIAL = LinearSegmentedColormap.from_list("house_blue", ["#ffffff", BLUE, DARK_BLUE])
DIVERGING = LinearSegmentedColormap.from_list("house_diverging", [ORANGE, "#ffffff", BLUE])

LETTER_SIZE = 20
DPI = 600


def apply() -> None:
    """Set the shared rcParams. Called on import, so importing is enough."""
    matplotlib.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,          # editable text in vector output
        "font.size": 12,
        "axes.labelsize": 14, "axes.titlesize": 13,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "legend.fontsize": 11, "legend.frameon": False,
        "axes.grid": False,
        "axes.prop_cycle": matplotlib.cycler(color=list(SERIES)),
        "lines.linewidth": 2.0, "lines.markersize": 5,
        "image.cmap": "house_blue",
    })


def panel_letter(ax, letter, x=-0.16, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=LETTER_SIZE, fontweight="bold",
            va="top", ha="left")


def save(fig, path_without_suffix):
    """Write ``<stem>.png`` (600 dpi) and ``<stem>.pdf``; return the PNG path."""
    for suffix in ("png", "pdf"):
        fig.savefig("%s.%s" % (path_without_suffix, suffix), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return "%s.png" % path_without_suffix


matplotlib.colormaps.register(SEQUENTIAL, force=True)
matplotlib.colormaps.register(DIVERGING, force=True)
apply()
