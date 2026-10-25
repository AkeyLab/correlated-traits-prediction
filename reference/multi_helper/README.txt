Equation 15 supplementary figure, version 2

PANEL A
-------
Input:
  eq15_validation_simulation.csv

Contains the original 300 validation simulations:
100 architectures each for m = 2, 5, and 10 helpers, N = 50,000 each.

The R plotting script fits one linear regression across all 300 architectures:
  Simulated R^2_gain = 0.0078 + 0.9861 * Predicted R^2_gain
  R^2 = 0.99701

PANEL B
-------
Inputs:
  eq15_panelB_simulation_replicates.csv
  eq15_panelB_simulation_summary.csv

Fixed parameters:
  m = 5 helpers
  target h^2 = 0.30
  helper h^2 = 0.20
  baseline R^2 = 0.10
  target-helper genetic correlation = 0
  target-helper environmental correlation = 0.25

The common helper-helper genetic and environmental correlation is varied
from 0 to 0.95 in increments of 0.05.

For every correlation value:
  100 replicate datasets were simulated
  N = 20,000 individuals per replicate

The plotted 95% simulation interval is the empirical 2.5th-97.5th percentile
range across the 100 replicate simulations. It is NOT a confidence interval
for the simulation mean.

The analytical line is Equation 15 evaluated from the generating parameters.

R PLOTTING
----------
plot_eq15_supplementary_figure_v2.R

Panel A:
  points = simulations
  dashed line = y=x
  solid black line = regression across all simulations
  annotation = fitted equation + R^2

Panel B:
  ribbon = empirical 95% simulation interval
  points = mean simulated R^2_gain
  solid line = analytical Equation 15 prediction

PLOTTING IN PYTHON
------------------
scripts/figureS_multi_helper.py draws the same figure in the house style
(used in the paper). The R script above uses the same colours and labels.
