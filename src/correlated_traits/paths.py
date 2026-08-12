"""Output locations.

Everything is written under ``results/`` at the repository root. Set the
``CORRELATED_TRAITS_RESULTS`` environment variable to redirect output elsewhere
(useful on a cluster where the repository lives on a read-only filesystem).

No absolute paths are hard-coded anywhere in this repository.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["REPO_ROOT", "RESULTS", "FIGURES", "ARRAYS", "TABLES", "ensure_output_dirs"]

REPO_ROOT = Path(__file__).resolve().parents[2]

RESULTS = Path(os.environ.get("CORRELATED_TRAITS_RESULTS", REPO_ROOT / "results"))
FIGURES = RESULTS / "figures"
ARRAYS = RESULTS / "arrays"
TABLES = RESULTS / "tables"


def ensure_output_dirs() -> None:
    """Create the output directories if they do not already exist."""
    for directory in (FIGURES, ARRAYS, TABLES):
        directory.mkdir(parents=True, exist_ok=True)
