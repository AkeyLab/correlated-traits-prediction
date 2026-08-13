"""Theory and simulation code for the correlated-traits prediction framework.

This package exists to reproduce the manuscript's figures. That is why the
per-figure modules, the ``correlated-traits`` command and the published
reference arrays all live inside it rather than beside it: one install, from
GitHub and without a clone, gives you everything needed to regenerate the
results and check them against what was published.

Package layout
    ``theory``    closed-form expressions for the relative predictive gain
    ``simulate``  the two-trait simulation model and the heatmap panel driver
    ``plotting``  shared figure style
    ``paths``     where outputs are written, and the packaged reference arrays
    ``figures``   one module per published figure, each exposing ``main()``
    ``cli``       the command that drives them
    ``selftest``  the checks behind ``correlated-traits test``

The manuscript's parameter grids are exported here because four figures and
three verification scripts all have to agree on them.
"""

from __future__ import annotations

import numpy as np

from . import paths, plotting, simulate, theory
from .simulate import SimulationConfig

__version__ = "1.0.0"

#: Correlation axis of every heatmap: 21 values from -1 to 1.
RHO_GRID = np.linspace(-1, 1, 21)

#: Helper-heritability axis of every heatmap: 9 values from 0.1 to 0.9.
H2_HELPER_GRID = np.round(np.linspace(0.1, 0.9, 9), 2)

#: Seed used for the published heatmap figures (3, 4, S2, S3).
HEATMAP_SEED = 20260717

#: Seed used for the published theory-validation figure (2).
VALIDATION_SEED = 42

#: Seed used for the independent replications behind Figure S1. Chosen to differ
#: from HEATMAP_SEED on purpose: these are independent replications, not
#: bit-for-bit reproductions.
REPLICATION_SEED = 20260803

__all__ = [
    "RHO_GRID", "H2_HELPER_GRID",
    "HEATMAP_SEED", "VALIDATION_SEED", "REPLICATION_SEED",
    "SimulationConfig", "paths", "plotting", "simulate", "theory", "__version__",
]
