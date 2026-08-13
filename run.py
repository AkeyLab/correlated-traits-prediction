#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["correlated-traits-prediction"]
#
# [tool.uv.sources]
# correlated-traits-prediction = { path = ".", editable = true }
# ///
"""Zero-install entry point for a bare clone.

    uv run ./run.py quick

uv reads the metadata block above, installs this repository into a throwaway
environment and dispatches to the same CLI as the ``correlated-traits`` console
script. If you have already run ``uv sync`` or ``pip install .``, use that
console script instead -- this file is only here so the first run needs nothing.
"""

from correlated_traits.cli import main

raise SystemExit(main())
