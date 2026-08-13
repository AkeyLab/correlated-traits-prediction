"""Pytest wrapper over :mod:`correlated_traits.selftest`.

The substance -- the verbatim transcription of the original implementation and
the bit-for-bit assertions -- lives in the package so that ``correlated-traits
test`` can run it without pytest and without a source checkout. This file exists
so ``pytest tests/`` picks the same checks up.

Run with:  pytest tests/ -v      (or: correlated-traits test)
"""

from __future__ import annotations

from correlated_traits.selftest import (
    check_closed_form_matches_least_squares,
    check_package_matches_original,
)


def test_package_matches_original_bit_for_bit():
    check_package_matches_original()


def test_closed_form_matches_simulation_identity():
    check_closed_form_matches_least_squares()
