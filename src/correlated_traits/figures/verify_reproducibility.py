"""Compare regenerated simulation arrays against the published ones.

The published arrays behind the figures in the paper ship inside the package,
under ``correlated_traits/reference/arrays/``. This module compares whatever is
currently in ``results/arrays/`` against them and reports, per array, whether
the match is bit-for-bit.

A bit-for-bit match is the expected outcome on the same NumPy version, because
the simulations are seeded. If the values agree to within floating-point noise
but are not identical, that almost always means a different BLAS or NumPy build
reordered a reduction; the science is unaffected and this says so rather than
failing. A large disagreement means something real has changed.

Exit status is 0 when every regenerated array matches to at least floating-point
tolerance, and 1 otherwise, so this is usable in continuous integration.

Usage:
    correlated-traits verify
"""

from __future__ import annotations

import numpy as np

from correlated_traits import paths

#: Arrays behind each published figure, and the command that regenerates them.
EXPECTED = {
    "Fig3_scenario1_highh_A.npy": "correlated-traits figure3-4",
    "Fig3_scenario1_highh_B.npy": "correlated-traits figure3-4",
    "Fig4_scenario2_lowh_A.npy": "correlated-traits figure3-4",
    "Fig4_scenario2_lowh_B.npy": "correlated-traits figure3-4",
    "S2_scenario3_lowh_goodbase_A.npy": "correlated-traits figureS2-S3",
    "S2_scenario3_lowh_goodbase_B.npy": "correlated-traits figureS2-S3",
    "S3_scenario4_highh_goodbase_A.npy": "correlated-traits figureS2-S3",
    "S3_scenario4_highh_goodbase_B.npy": "correlated-traits figureS2-S3",
}

TOLERANCE = 1e-10


def main() -> int:
    print("Comparing %s\n     against %s\n" % (paths.ARRAYS, paths.reference_location()))
    print("%-38s %-12s %12s" % ("array", "status", "max |diff|"))
    print("-" * 64)

    identical = close = missing = failed = no_reference = 0

    for name, producer in sorted(EXPECTED.items()):
        regenerated_path = paths.ARRAYS / name

        if not paths.has_reference_array(name):
            print("%-38s %-12s %12s" % (name, "no ref", "-"))
            no_reference += 1
            continue
        if not regenerated_path.exists():
            print("%-38s %-12s %12s   run %s" % (name, "not run", "-", producer))
            missing += 1
            continue

        reference = paths.reference_array(name)
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
          "%d reference missing"
          % (identical, close, missing, failed, no_reference))

    # Checked before anything else: with no references to compare against, every
    # other counter stays at zero, which would otherwise read as a clean pass.
    if no_reference:
        print("\nFAIL: %d reference array(s) are missing. The published arrays ship inside "
              "the package, so this usually means a broken or incomplete install."
              % no_reference)
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
