# Reproduce the figures. `make figures` builds everything in dependency order.
#
# PYTHON can be overridden:  make figures PYTHON=/path/to/python

PYTHON ?= python
SCRIPTS = scripts
RESULTS ?= results

.PHONY: help quick figures heatmaps replications supplementary restyle test clean

help:
	@echo "make quick        The test suite                       (~10 seconds)"
	@echo "make figures      Every figure, in dependency order    (~1 hour)"
	@echo "make heatmaps     Figures 3, 4, S1, S2 only            (~50 minutes)"
	@echo "make replications Uncertainty estimates for Figure 2D-F (~40 minutes)"
	@echo "make supplementary Supplementary Figures 3 and 4       (~3 minutes)"
	@echo "make restyle      Redraw Figures 3, 4, S1, S2 from arrays (seconds)"
	@echo "make test         Bit-for-bit regression test          (~10 seconds)"
	@echo "make verify       Compare output against reference/    (~1 second)"
	@echo "make clean        Delete everything under $(RESULTS)/"

quick: test

test:
	$(PYTHON) tests/test_equivalence.py

heatmaps:
	$(PYTHON) $(SCRIPTS)/figure3_figure4_heatmaps.py
	$(PYTHON) $(SCRIPTS)/figureS1_figureS2_heatmaps.py

# Supplementary Figures 3 and 4 (S3 needs the heatmap arrays; S4 simulates, ~3 minutes)
supplementary:
	$(PYTHON) $(SCRIPTS)/figureS_phenotypic_correlation.py
	$(PYTHON) $(SCRIPTS)/figureS_multi_helper.py

# Redraw Figures 3, 4, S1, S2 from saved arrays without re-simulating (seconds)
restyle:
	$(PYTHON) $(SCRIPTS)/render_heatmaps.py

replications:
	$(PYTHON) $(SCRIPTS)/replicate_figure3A_column.py
	$(PYTHON) $(SCRIPTS)/replicate_figure3A_figure4A_panels.py

# Figure 2 panels D-F need the simulated values (heatmaps) and the standard
# errors (replications), so Figure 2 must come last.
figures: test
	$(MAKE) heatmaps PYTHON=$(PYTHON)
	$(MAKE) replications PYTHON=$(PYTHON)
	$(PYTHON) $(SCRIPTS)/figure2_theory_validation.py
	$(MAKE) supplementary PYTHON=$(PYTHON)
	$(PYTHON) $(SCRIPTS)/verify_reproducibility.py

verify:
	$(PYTHON) $(SCRIPTS)/verify_reproducibility.py

clean:
	find $(RESULTS) -type f ! -name '.gitkeep' -delete
