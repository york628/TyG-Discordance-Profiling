# =============================================================================
# 02_biomarker_profiles.R
# TyG-Biomarker Discordance Profiling
# Step 2: Cluster weight bar charts + biomarker residual scatter plots
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
# Output  : TyG_cluster_overview.png / .pdf
# =============================================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)

OUTPUT_DIR <- "."

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))

# =============================================================================
# Constants
# =============================================================================

BIOMARKERS    <- c("bmi", "sbp", "dbp", "egfr", "crp", "hdl", "ldl", "chol")
SEX_GROUPS    <- c("Female", "Male")
CLUSTER_ORDER <- c("BC", "DHT", "DRI", "DHL", "DOB", "DPL", "DIS")

BIOMARKER_ORDER  <- rev(c("BMI", "SBP", "DBP", "EGFR", "CRP", "HDL", "LDL", "CHOL"))
BIOMARKER_LABELS <- c(
  BMI = "BMI", SBP = "SBP", DBP = "DBP", EGFR = "eGFR",
  CRP = "CRP", HDL = "HDL", LDL = "LDL", CHOL = "CHOL"
)

CLUSTER_COLORS <- c(
  BC  = "#3d7d8a", DHT = "#5ea898", DRI = "#84c49a",
  DHL = "#b8d990", DOB = "#f0d882", DPL = "#eba85c", DIS = "#e08a60"
)

# Mapping: semantic label → raw column name in clustering results
CLUSTER_MAP <- list(
  Female = list(
    BC = "cluster_0",  DHT = "cluster_3",  DRI = "cluster_2",
    DHL = "cluster_7", DOB = "cluster_12", DPL = "cluster_1", DIS = "cluster_9"
  ),
  Male = list(
    BC = "cluster_0",  DHT = "cluster_3",  DRI = "cluster_2",
    DHL = "cluster_6", DOB = "cluster_12", DPL = "cluster_5", DIS = "cluster_13"
  )
)

# =============================================================================
# Data extraction
# =============================================================================

extract_weights <- function(sex_name) {
  clusters <- complete_results$clustering_results[[sex_name]]$clusters
  mapping  <- CLUSTER_MAP[[sex_name]]
  map_dfr(names(mapping), function(label) {
    cl_name <- mapping[[label]]
    if (!cl_name %in% names(clusters)) return(NULL)
    data.frame(
      sex     = sex_name,
      profile = label,
      weight  = clusters[[cl_name]]$weight * 100
    )
  })
}

extract_centers <- function(sex_name) {
  clusters <- complete_results$clustering_results[[sex_name]]$clusters
  mapping  <- CLUSTER_MAP[[sex_name]]
  map_dfr(names(mapping), function(label) {
    cl_name <- mapping[[label]]
    if (!cl_name %in% names(clusters)) return(NULL)
    center <- clusters[[cl_name]]$center
    map_dfr(BIOMARKERS, function(bm) {
      if (!bm %in% names(center)) return(NULL)
      data.frame(
        sex       = sex_name,
        profile   = label,
        biomarker = toupper(bm),
        value     = center[[bm]]
      )
    })
  })
}

weight_data <- bind_rows(map(SEX_GROUPS, extract_weights)) |>
  mutate(profile = factor(profile, levels = CLUSTER_ORDER))

center_data <- bind_rows(map(SEX_GROUPS, extract_centers)) |>
  mutate(
    biomarker = factor(biomarker, levels = BIOMARKER_ORDER),
    profile   = factor(profile,   levels = CLUSTER_ORDER)
  )

# =============================================================================
# Bar chart: cluster proportions
# =============================================================================

make_bar <- function(sex_name, show_ylab = TRUE) {
  weight_data |>
    filter(sex == sex_name) |>
    ggplot(aes(x = profile, y = weight, fill = profile)) +
    geom_col(width = 0.55) +
    scale_fill_manual(values = CLUSTER_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = sex_name,
      x     = "Profile",
      y     = if (show_ylab) "Percentage (%)" else NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title   = element_text(size = 9, hjust = 0.5, face = "bold"),
      axis.text.x  = element_text(size = 8, face = "bold"),
      axis.text.y  = element_text(size = 7),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "grey70", linewidth = 0.4),
      aspect.ratio = 0.65,
      plot.margin  = margin(3, 4, 2, 4, "mm")
    )
}

# =============================================================================
# Scatter plot: cluster centroid residuals
# =============================================================================

scatter_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title         = element_blank(),
    axis.text.x        = element_text(size = 7.5, face = "bold", color = "black"),
    axis.text.y        = element_text(size = 6.5, color = "black"),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.text.y.right = element_text(
      size = 8.5, face = "bold", angle = 0, hjust = 0.5,
      margin = margin(l = 3, r = 3)
    ),
    strip.background = element_rect(fill = "grey90", color = "grey70", linewidth = 0.3),
    strip.placement  = "outside",
    panel.spacing.y  = unit(0.03, "lines"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(2, 2, 2, 2, "mm")
  )

make_scatter <- function(sex_name, show_strip = TRUE) {
  p <- center_data |>
    filter(sex == sex_name) |>
    ggplot(aes(x = profile, y = value, color = profile)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = CLUSTER_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    facet_grid(
      biomarker ~ ., scales = "free_y",
      labeller = labeller(biomarker = BIOMARKER_LABELS)
    ) +
    labs(title = sex_name) +
    scatter_theme +
    theme(
      plot.title = element_text(
        size = 10, face = "bold", hjust = 0.5, margin = margin(b = 3)
      )
    )
  if (!show_strip) p <- p + theme(strip.text.y = element_blank())
  p
}

# =============================================================================
# Assemble and save
# =============================================================================

bar_row <- plot_grid(
  make_bar("Female", TRUE), make_bar("Male", FALSE),
  nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb"
)

scat_pair <- plot_grid(
  make_scatter("Female", TRUE), make_scatter("Male", FALSE),
  nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb"
)

y_label  <- ggdraw() +
  draw_label("Residuals \u2013 s.d. units", size = 8, angle = 90)
scat_row <- plot_grid(y_label, scat_pair, nrow = 1, rel_widths = c(0.03, 1))

bar_row_aligned <- plot_grid(ggdraw(), bar_row, nrow = 1, rel_widths = c(0.03, 1))
final_plot      <- plot_grid(bar_row_aligned, scat_row, ncol = 1, rel_heights = c(1.8, 4.5))

ggsave(file.path(OUTPUT_DIR, "TyG_cluster_overview.png"),
       final_plot, width = 8, height = 10, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "TyG_cluster_overview.pdf"),
       final_plot, width = 8, height = 10, bg = "white")

cat("Done. Saved: TyG_cluster_overview.png / .pdf\n")
