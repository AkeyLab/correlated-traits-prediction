"""Output locations and packaged reference data.

Results are written under ``results/`` in the current working directory, so
running the CLI from a checkout puts them where they have always gone. Three
ways to redirect that, highest precedence first:

    ``correlated-traits --results-dir PATH``   per invocation
    ``set_results_root(path)``                 from Python
    ``CORRELATED_TRAITS_RESULTS``              environment variable

The paths resolve lazily, at the moment they are read, so setting the
environment variable after importing this module still works.

The published reference arrays ship inside the package rather than sitting at
the repository root, so ``correlated-traits verify`` works from an installed
copy and from any directory.

No absolute paths are hard-coded anywhere in this repository.
"""

from __future__ import annotations

import io
import os
from importlib import resources
from pathlib import Path

import numpy as np

__all__ = [
    "RESULTS", "FIGURES", "ARRAYS", "TABLES",
    "set_results_root", "ensure_output_dirs",
    "reference_array", "has_reference_array", "reference_location",
]

_REFERENCE_PACKAGE = "correlated_traits"
_REFERENCE_SUBPATH = ("reference", "arrays")

_override: Path | None = None


def set_results_root(path) -> None:
    """Point every output directory at ``path`` for the rest of the process."""
    global _override
    _override = Path(path).expanduser().resolve()


def _results_root() -> Path:
    if _override is not None:
        return _override
    environment = os.environ.get("CORRELATED_TRAITS_RESULTS")
    if environment:
        return Path(environment).expanduser()
    return Path.cwd() / "results"


_SUBDIRECTORIES = {"FIGURES": "figures", "ARRAYS": "arrays", "TABLES": "tables"}


def __getattr__(name: str) -> Path:
    """Resolve RESULTS/FIGURES/ARRAYS/TABLES on access rather than on import."""
    if name == "RESULTS":
        return _results_root()
    if name in _SUBDIRECTORIES:
        return _results_root() / _SUBDIRECTORIES[name]
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(list(globals()) + ["RESULTS"] + list(_SUBDIRECTORIES))


def ensure_output_dirs() -> None:
    """Create the output directories if they do not already exist."""
    root = _results_root()
    for subdirectory in _SUBDIRECTORIES.values():
        (root / subdirectory).mkdir(parents=True, exist_ok=True)


def _reference_traversable(name: str):
    directory = resources.files(_REFERENCE_PACKAGE)
    for part in _REFERENCE_SUBPATH:
        directory = directory / part
    return directory / name


def has_reference_array(name: str) -> bool:
    """Whether ``name`` is among the reference arrays shipped in the package."""
    try:
        return _reference_traversable(name).is_file()
    except (FileNotFoundError, ModuleNotFoundError):
        return False


def reference_array(name: str) -> np.ndarray:
    """Load a published reference array from the package's own data."""
    return np.load(io.BytesIO(_reference_traversable(name).read_bytes()))


def reference_location() -> str:
    """Human-readable location of the packaged reference arrays, for reporting."""
    directory = resources.files(_REFERENCE_PACKAGE)
    for part in _REFERENCE_SUBPATH:
        directory = directory / part
    return str(directory)
