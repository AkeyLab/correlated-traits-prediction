"""Command-line entry point: ``correlated-traits <command>``.

Every command regenerates part of the manuscript and writes under ``results/``.
The composite commands exist because the figures have prerequisites -- see
``figures`` below, where the ordering is load-bearing.

Figure modules are imported inside the handlers, not at module level. Several of
them pull in matplotlib through :mod:`correlated_traits.plotting`, which forces
the Agg backend at import time, and ``--help`` should not pay for that.
"""

from __future__ import annotations

import argparse
import shutil
import sys

__all__ = ["main"]


def _run_module(module_name: str):
    """Import ``correlated_traits.figures.<module_name>`` and call its ``main()``."""
    from importlib import import_module

    module = import_module("correlated_traits.figures." + module_name)
    return module.main()


def _selftest() -> None:
    from correlated_traits import selftest

    selftest.run()


# Single-step commands: name -> (help text, callable).
STEPS = {
    "test": ("Bit-for-bit regression check (~10 seconds)",
             _selftest),
    "figure2": ("Figure 2, theory validation (seconds)",
                lambda: _run_module("figure2_theory_validation")),
    "figure3-4": ("Figures 3 and 4, poor-baseline heatmaps (30-60 minutes)",
                  lambda: _run_module("figure3_figure4_heatmaps")),
    "figureS2-S3": ("Figures S2 and S3, strong-baseline heatmaps (30-60 minutes)",
                    lambda: _run_module("figureS2_figureS3_heatmaps")),
    "figureS1": ("Figure S1, closed form vs simulation (seconds; needs figure3-4 "
                 "and replicate-panel first)",
                 lambda: _run_module("figureS1_theory_vs_simulation")),
    "replicate-column": ("Replicate one column of Figure 3A (a few minutes)",
                         lambda: _run_module("replicate_figure3A_column")),
    "replicate-panel": ("Replicate all of Figure 3A with error bars (10-30 minutes)",
                        lambda: _run_module("replicate_figure3A_panel")),
    "verify": ("Compare results/arrays against the published arrays (~1 second)",
               lambda: _run_module("verify_reproducibility")),
}

# Composite commands: name -> (help text, ordered list of step names).
#
# The order in `figures` matters. Figure S1 consumes
# results/arrays/Fig3_scenario1_highh_A.npy (from figure3-4) and
# results/tables/replication_fig3A_panel.npz (from replicate-panel), so it must
# come after both.
PIPELINES = {
    "quick": ("Figure 2 and the regression check (~30 seconds)",
              ["test", "figure2"]),
    "heatmaps": ("Figures 3, 4, S2 and S3 (1-2 hours)",
                 ["figure3-4", "figureS2-S3"]),
    "replications": ("Uncertainty estimates behind Figure S1 (15-35 minutes)",
                     ["replicate-column", "replicate-panel"]),
    "figures": ("Every figure, in dependency order (1.5-2.5 hours)",
                ["test", "figure2", "figure3-4", "figureS2-S3",
                 "replicate-column", "replicate-panel", "figureS1", "verify"]),
}


def _version() -> str:
    """Package version, without importing the package (and so without matplotlib)."""
    from importlib import metadata

    try:
        return metadata.version("correlated-traits-prediction")
    except metadata.PackageNotFoundError:
        from correlated_traits import __version__

        return __version__


def _clean(results_root, assume_yes: bool) -> int:
    if not results_root.exists():
        print("Nothing to clean: %s does not exist." % results_root)
        return 0
    if not assume_yes:
        if not sys.stdin.isatty():
            print("Refusing to delete %s without --yes." % results_root, file=sys.stderr)
            return 1
        reply = input("Delete everything under %s? [y/N] " % results_root)
        if reply.strip().lower() not in ("y", "yes"):
            print("Cancelled.")
            return 1
    shutil.rmtree(results_root, ignore_errors=True)
    print("Deleted %s" % results_root)
    return 0


def build_parser() -> argparse.ArgumentParser:
    # --results-dir is accepted on either side of the subcommand, because both
    # readings are natural. SUPPRESS is what makes that safe: without it the
    # subparser's default would overwrite a value given before the subcommand.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--results-dir", metavar="PATH", default=argparse.SUPPRESS,
        help="Write results here instead of ./results.",
    )

    parser = argparse.ArgumentParser(
        prog="correlated-traits",
        description="Regenerate the figures for the correlated-traits prediction paper.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[common],
        epilog="Outputs go to ./results by default. Override with --results-dir or "
               "the CORRELATED_TRAITS_RESULTS environment variable.",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + _version())
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    for name, (help_text, steps) in PIPELINES.items():
        subparsers.add_parser(name, help=help_text, parents=[common],
                              description="%s\n\nRuns: %s" % (help_text, ", ".join(steps)),
                              formatter_class=argparse.RawDescriptionHelpFormatter)
    for name, (help_text, _) in STEPS.items():
        subparsers.add_parser(name, help=help_text, description=help_text, parents=[common])

    clean = subparsers.add_parser("clean", help="Delete everything under results/",
                                  parents=[common])
    clean.add_argument("--yes", action="store_true",
                       help="Do not ask for confirmation.")

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    # Must happen before any figure module is imported, so that paths.ARRAYS and
    # friends resolve against the requested root.
    from correlated_traits import paths

    results_dir = getattr(args, "results_dir", None)
    if results_dir:
        paths.set_results_root(results_dir)

    if args.command == "clean":
        return _clean(paths.RESULTS, args.yes)

    steps = PIPELINES[args.command][1] if args.command in PIPELINES else [args.command]

    for step in steps:
        if len(steps) > 1:
            print("\n>>> correlated-traits %s" % step, flush=True)
        status = STEPS[step][1]()
        if status:
            return int(status)
    return 0
