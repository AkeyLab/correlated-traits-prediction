#!/usr/bin/env python3
"""Compare regenerated simulation arrays against the published ones.

``reference/arrays/`` holds the arrays behind the figures in the paper. This
script compares whatever is currently in ``results/arrays/`` against them and
reports, per array, whether the match is bit-for-bit.

A bit-for-bit match is the expected outcome on the same NumPy version, because
the simulations are seeded. If the values agree to within floating-point noise
but are not identical, that almost always means a different BLAS or NumPy build
reordered a reduction; the science is unaffected and the script says so rather
than failing. A large disagreement means something real has changed.

Exit status is 0 when every regenerated array matches to at least floating-point
tolerance, and 1 otherwise, so this is usable in continuous integration.

Usage:
    python scripts/verify_reproducibility.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from correlated_traits import paths  # noqa: E402

REFERENCE_DIR = paths.REPO_ROOT / "reference" / "arrays"

#: Arrays behind each published figure, and the script that regenerates them.
EXPECTED = {
    "Fig3_scenario1_highh_A.npy": "figure3_figure4_heatmaps.py",
    "Fig3_scenario1_highh_B.npy": "figure3_figure4_heatmaps.py",
    "Fig4_scenario2_lowh_A.npy": "figure3_figure4_heatmaps.py",
    "Fig4_scenario2_lowh_B.npy": "figure3_figure4_heatmaps.py",
    "S2_scenario3_lowh_goodbase_A.npy": "figureS2_figureS3_heatmaps.py",
    "S2_scenario3_lowh_goodbase_B.npy": "figureS2_figureS3_heatmaps.py",
    "S3_scenario4_highh_goodbase_A.npy": "figureS2_figureS3_heatmaps.py",
    "S3_scenario4_highh_goodbase_B.npy": "figureS2_figureS3_heatmaps.py",
}

TOLERANCE = 1e-10


def main() -> int:
    if not REFERENCE_DIR.is_dir():
        print("No reference directory at %s" % REFERENCE_DIR)
        return 1

    print("Comparing %s\n     against %s\n" % (paths.ARRAYS, REFERENCE_DIR))
    print("%-38s %-12s %12s" % ("array", "status", "max |diff|"))
    print("-" * 64)

    identical = close = missing = failed = no_reference = 0

    for name, producer in sorted(EXPECTED.items()):
        reference_path = REFERENCE_DIR / name
        regenerated_path = paths.ARRAYS / name

        if not reference_path.exists():
            print("%-38s %-12s %12s" % (name, "NO REF", "-"))
            no_reference += 1
            continue
        if not regenerated_path.exists():
            print("%-38s %-12s %12s   run %s" % (name, "not run", "-", producer))
            missing += 1
            continue

        reference = np.load(reference_path)
        regenerated = np.load(regenerated_path)

        if reference.shape != regenerated.shape:
            print("%-38s %-12s %12s   %s vs %s"
                  % (name, "SHAPE", "-", reference.shape, regenerated.shape))
            failed += 1
            continue

        difference = float(np.abs(reference - regenerated).max())
        if np.array_equal(reference, regenerated):
            print("%-38s %-12s %12s" % (name, "identical", "0"))
            identical += 1
        elif difference < TOLERANCE:
            print("%-38s %-12s %12.3e" % (name, "close", difference))
            close += 1
        else:
            print("%-38s %-12s %12.3e" % (name, "DIFFERS", difference))
            failed += 1

    print("-" * 64)
    print("%d identical, %d within tolerance, %d not yet regenerated, %d differing, "
          "%d with no reference"
          % (identical, close, missing, failed, no_reference))

    if no_reference:
        print("\nFAIL: %d of %d reference arrays are missing from %s, so there was nothing "
              "to compare against. An empty reference directory must not report success."
              % (no_reference, len(EXPECTED), REFERENCE_DIR))
        return 1
    if not (identical or close):
        print("\nFAIL: no array was actually compared.")
        return 1
    if failed:
        print("\nFAIL: some arrays differ by more than %g." % TOLERANCE)
        return 1
    if missing:
        print("\nIncomplete: regenerate the missing arrays, then re-run.")
        return 1
    if close:
        print("\nPASS with floating-point differences, most likely a different "
              "NumPy or BLAS build. Conclusions are unaffected.")
        return 0
    print("\nPASS: every array reproduces bit-for-bit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
