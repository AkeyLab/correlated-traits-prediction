# Reproduce the figures. `make figures` builds everything in dependency order.
#
# PYTHON can be overridden:  make figures PYTHON=/path/to/python

PYTHON ?= python
SCRIPTS = scripts
RESULTS ?= results

.PHONY: help quick figures heatmaps replications test clean

help:
	@echo "make quick        Figure 2 and the test suite          (~30 seconds)"
	@echo "make figures      Every figure, in dependency order    (~1 hour)"
	@echo "make heatmaps     Figures 3, 4, S2, S3 only            (~50 minutes)"
	@echo "make replications Uncertainty estimates for Figure S1  (~22 minutes)"
	@echo "make test         Bit-for-bit regression test          (~10 seconds)"
	@echo "make verify       Compare output against reference/    (~1 second)"
	@echo "make clean        Delete everything under $(RESULTS)/"

quick: test
	$(PYTHON) $(SCRIPTS)/figure2_theory_validation.py

test:
	$(PYTHON) tests/test_equivalence.py

heatmaps:
	$(PYTHON) $(SCRIPTS)/figure3_figure4_heatmaps.py
	$(PYTHON) $(SCRIPTS)/figureS2_figureS3_heatmaps.py

replications:
	$(PYTHON) $(SCRIPTS)/replicate_figure3A_column.py
	$(PYTHON) $(SCRIPTS)/replicate_figure3A_panel.py

# Figure S1 needs the simulated values (heatmaps) and the standard errors
# (replications), so it must come last.
figures: test
	$(PYTHON) $(SCRIPTS)/figure2_theory_validation.py
	$(MAKE) heatmaps PYTHON=$(PYTHON)
	$(MAKE) replications PYTHON=$(PYTHON)
	$(PYTHON) $(SCRIPTS)/figureS1_theory_vs_simulation.py
	$(PYTHON) $(SCRIPTS)/verify_reproducibility.py

verify:
	$(PYTHON) $(SCRIPTS)/verify_reproducibility.py

clean:
	find $(RESULTS) -type f ! -name '.gitkeep' -delete
