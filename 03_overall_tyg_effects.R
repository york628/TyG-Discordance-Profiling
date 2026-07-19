# =============================================================================
# 03_overall_tyg_effects.R
# TyG-Biomarker Discordance Profiling
# Step 3: Forest plots of overall TyG → biomarker effects (two-row layout)
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
# Output  : TyG_overall_effects.png / .pdf
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

OUTPUT_DIR <- "."

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))

# =============================================================================
# Constants
# =============================================================================

BIOMARKER_LABELS <- c(
  bmi  = "BMI\n(kg/m\u00b2)",
  sbp  = "SBP\n(mmHg)",
  dbp  = "DBP\n(mmHg)",
  egfr = "eGFR\n(mL/min/1.73m\u00b2)",
  crp  = "CRP\n(mg/L)",
  hdl  = "HDL\n(mmol/L)",
  ldl  = "LDL\n(mmol/L)",
  chol = "CHOL\n(mmol/L)"
)

ROW1 <- c("bmi", "sbp", "dbp", "egfr")
ROW2 <- c("crp", "hdl", "ldl", "chol")

# =============================================================================
# Data preparation
# =============================================================================

plot_data <- complete_results$bmi_coefs |>
  filter(term == "tyg") |>
  mutate(
    is_sig            = !(lowerCI < 0 & upperCI > 0),
    sex_label         = factor(sex, levels = c("Female", "Male")),
    biomarker_display = factor(
      Biomarker,
      levels = c(ROW1, ROW2),
      labels = BIOMARKER_LABELS[c(ROW1, ROW2)]
    ),
    row_group = ifelse(Biomarker %in% ROW1, "Row1", "Row2")
  ) |>
  filter(!is.na(biomarker_display))

# =============================================================================
# Forest plot function
# =============================================================================

make_forest <- function(df, show_title = FALSE) {
  ggplot(df, aes(x = Estimate, y = sex_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_linerange(aes(xmin = lowerCI, xmax = upperCI),
                   linewidth = 1.0, color = "black") +
    geom_point(aes(fill = is_sig), shape = 21, size = 3.5,
               stroke = 0.8, color = "black") +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    facet_grid(. ~ biomarker_display, scales = "free_x", space = "free_y") +
    labs(
      x     = if (show_title) "Estimate (95% CI)" else NULL,
      y     = NULL,
      title = if (show_title) "Overall TyG effects on biomarkers" else NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x        = element_text(size = 7),
      axis.text.y        = element_text(size = 9, face = "bold"),
      axis.title.x       = element_text(size = 9, margin = margin(t = 6)),
      strip.text.x       = element_text(size = 8, face = "bold"),
      strip.background   = element_rect(fill = "grey85", color = "black", linewidth = 0.4),
      plot.title         = element_text(size = 11, hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border       = element_rect(color = "black", linewidth = 0.4),
      panel.spacing.x    = unit(0.8, "lines"),
      plot.margin        = margin(4, 6, 4, 6)
    )
}

# =============================================================================
# Assemble and save
# =============================================================================

p_row1 <- make_forest(
  filter(plot_data, row_group == "Row1") |> mutate(biomarker_display = droplevels(biomarker_display)),
  show_title = TRUE
)
p_row2 <- make_forest(
  filter(plot_data, row_group == "Row2") |> mutate(biomarker_display = droplevels(biomarker_display)),
  show_title = FALSE
)

final_plot <- p_row1 / p_row2 + plot_layout(heights = c(1, 1))

ggsave(file.path(OUTPUT_DIR, "TyG_overall_effects.png"),
       final_plot, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "TyG_overall_effects.pdf"),
       final_plot, width = 10, height = 6, bg = "white")

cat("Done. Saved: TyG_overall_effects.png / .pdf\n")
