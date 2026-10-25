# Supplementary Figure: validation and illustration of Equation 15
#
# Panel A:
#   Predicted vs. simulated R^2_gain across 300 randomly generated architectures
#   (100 each with 2, 5, and 10 helper traits).
#   A single least-squares regression is fit across all 300 simulations.
#
# Panel B:
#   Five-helper simulation in which only helper-helper correlation changes.
#   The solid line is the analytical prediction from Equation 15.
#   Points are mean simulated R^2_gain across 100 replicate datasets.
#   The ribbon is the empirical 95% simulation interval
#   (2.5th and 97.5th percentiles across replicates).
#
# Required packages:
# install.packages(c("ggplot2", "dplyr", "patchwork"))

library(ggplot2)
library(dplyr)
library(patchwork)
library(latex2exp)

# -------------------------------------------------------------------------
# Input files
# -------------------------------------------------------------------------

data_dir <- "."

validation <- read.csv(
  file.path(data_dir, "eq15_validation_simulation.csv")
)

panelB <- read.csv(
  file.path(data_dir, "eq15_panelB_simulation_summary.csv")
)

# -------------------------------------------------------------------------
# Aesthetic controls
# -------------------------------------------------------------------------
# Colours match the house palette (src/correlated_traits/style.py): blue, orange, green.

helper_colors <- c(
  "2 helpers"  = "#1f77b4",
  "5 helpers"  = "#ff7f0e",
  "10 helpers" = "#2ca02c"
)

helper_shapes <- c(
  "2 helpers"  = 16,
  "5 helpers"  = 17,
  "10 helpers" = 15
)

analytical_color <- "black"
simulation_color <- "#1f77b4"

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.tag = element_text(face = "bold", size = 13),
    legend.title = element_blank(),
    legend.position = "top"
  )

# =========================================================================
# PANEL A: Equation 15 predictions vs. simulated R^2_gain
# =========================================================================

validation <- validation %>%
  mutate(
    helper_set = factor(
      paste(m, "helpers"),
      levels = c("2 helpers", "5 helpers", "10 helpers")
    )
  )

# Fit one regression across all 300 simulated architectures.
fit_A <- lm(
  r2_gain_simulated ~ r2_gain_predicted,
  data = validation
)

fit_intercept <- unname(coef(fit_A)[1])
fit_slope <- unname(coef(fit_A)[2])
fit_r2 <- summary(fit_A)$r.squared

# Equation annotation.
# Change digits here if you want more/less precision.
fit_label <- sprintf(
  "y = %.3f + %.3fx\nR² = %.4f",
  fit_intercept,
  fit_slope,
  fit_r2
)

axis_max <- max(
  validation$r2_gain_predicted,
  validation$r2_gain_simulated
)

pA <- ggplot(
  validation,
  aes(
    x = r2_gain_predicted,
    y = r2_gain_simulated,
    color = helper_set,
    shape = helper_set
  )
) +

  # Expected relationship under perfect agreement.
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.6,
    linetype = "dashed",
    color = "grey50"
  ) +

  # Least-squares fit across all 300 architectures.
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +

  geom_point(
    size = 2.0,
    alpha = 0.70
  ) +

  annotate(
    "text",
    x = 0.05 * axis_max,
    y = 0.95 * axis_max,
    label = fit_label,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +

  scale_color_manual(values = helper_colors) +
  scale_shape_manual(values = helper_shapes) +

  coord_equal(
    xlim = c(0, axis_max * 1.03),
    ylim = c(0, axis_max * 1.03),
    expand = FALSE
  ) +

  labs(
    x = TeX("Predicted $R^2_{gain}$"),
    y = TeX("Simulated $R^2_{gain}$")
  ) +

  base_theme

# =========================================================================
# PANEL B: Helper-helper correlation
# =========================================================================
#
# At every x value:
#   - helper marginal relationships with the target are fixed
#   - only pairwise correlation among helpers changes
#   - the ribbon is the empirical 95% simulation interval
#   - points are the mean across simulation replicates
#   - the line is the Equation 15 analytical prediction

pB <- ggplot(
  panelB,
  aes(x = pairwise_helper_correlation)
) +

  # Empirical 95% simulation envelope.
  geom_ribbon(
    aes(
      ymin = r2_gain_sim_lower95,
      ymax = r2_gain_sim_upper95
    ),
    fill = simulation_color,
    alpha = 0.18
  ) +

  # Analytical prediction from Equation 15.
  geom_line(
    aes(
      y = r2_gain_analytical,
      linetype = "Equation 15"
    ),
    color = analytical_color,
    linewidth = 1.0
  ) +

  # Mean of the replicate simulations.
  geom_point(
    aes(
      y = r2_gain_sim_mean,
      shape = "Simulation mean"
    ),
    color = simulation_color,
    size = 2.1
  ) +

  scale_linetype_manual(
    values = c("Equation 15" = "solid"),
    name = NULL
  ) +

  scale_shape_manual(
    values = c("Simulation mean" = 16),
    name = NULL
  ) +

  labs(
    x = "Pairwise correlation among helper traits",
    y = TeX("$R^2_{gain}$")
  ) +

  base_theme +

  theme(
    legend.position = "top"
  )

# =========================================================================
# Combine panels
# =========================================================================

fig <- pA + pB +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 1))

print(fig)

# -------------------------------------------------------------------------
# Save figure
# -------------------------------------------------------------------------

ggsave(
  filename = "eq15_supplementary_figure_v2.pdf",
  plot = fig,
  width = 7.3,
  height = 3.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "eq15_supplementary_figure_v2.png",
  plot = fig,
  width = 7.3,
  height = 3.6,
  units = "in",
  dpi = 600
)

# Optional individual panels.
ggsave(
  "eq15_panel_A_v2.pdf",
  pA,
  width = 3.6,
  height = 3.5,
  units = "in"
)

ggsave(
  "eq15_panel_B_v2.pdf",
  pB,
  width = 3.6,
  height = 3.5,
  units = "in"
)

# -------------------------------------------------------------------------
# Print regression statistics for Panel A
# -------------------------------------------------------------------------

cat("\nPanel A regression:\n")
cat(sprintf(
  "Simulated R^2_gain = %.4f + %.4f * Predicted R^2_gain\n",
  fit_intercept,
  fit_slope
))
cat(sprintf("R^2 = %.5f\n", fit_r2))
