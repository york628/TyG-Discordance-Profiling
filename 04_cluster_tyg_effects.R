# =============================================================================
# 04_cluster_tyg_effects.R
# TyG-Biomarker Discordance Profiling
# Step 4: Cluster-specific TyG effects — weighted regression, forest plot,
#         and supplementary Excel table
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
# Output  : TyG_cluster_effects.png / .pdf
#           TyG_cluster_effects.csv
#           Supplementary_Table_Cluster_TyG_Effects.xlsx
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(openxlsx)

OUTPUT_DIR <- "."

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))

# =============================================================================
# Constants
# =============================================================================

BIOMARKERS    <- c("bmi", "sbp", "dbp", "egfr", "crp", "hdl", "ldl", "chol")
EXPOSURE      <- "tyg"
COVARIATES    <- c("age", "smoking")
SEX_GROUPS    <- c("Female", "Male")
CLUSTER_ORDER <- c("BC", "DHT", "DRI", "DHL", "DOB", "DPL", "DIS")

BIOMARKER_LABELS <- c(
  bmi  = "BMI\n(kg m\u207b\u00b2)",
  sbp  = "SBP\n(mmHg)",
  dbp  = "DBP\n(mmHg)",
  egfr = "eGFR\n(mL min\u207b\u00b9\n1.73m\u207b\u00b2)",
  crp  = "CRP\n(mg L\u207b\u00b9)",
  hdl  = "HDL\n(mmol L\u207b\u00b9)",
  ldl  = "LDL\n(mmol L\u207b\u00b9)",
  chol = "CHOL\n(mmol L\u207b\u00b9)"
)

BIOMARKER_UNITS <- c(
  bmi  = "kg/m\u00b2", sbp  = "mmHg",    dbp  = "mmHg",
  egfr = "mL/min/1.73m\u00b2",            crp  = "mg/L",
  hdl  = "mmol/L",    ldl  = "mmol/L",   chol = "mmol/L"
)

PHENOTYPE_LABELS <- c(
  BC  = "Baseline Concordant",
  DHT = "Discordant Hypertensive",
  DRI = "Discordant Renal Insufficiency",
  DHL = "Discordant High HDL-Lipid",
  DOB = "Discordant Obesity",
  DPL = "Discordant Pro-Atherogenic Lipid",
  DIS = "Discordant Inflammatory"
)

CLUSTER_MAP <- list(
  Female = list(
    BC  = "cluster_0", DHT = "cluster_3",  DRI = "cluster_2",
    DHL = "cluster_7", DOB = "cluster_12", DPL = "cluster_1", DIS = "cluster_9"
  ),
  Male = list(
    BC  = "cluster_0", DHT = "cluster_3",  DRI = "cluster_2",
    DHL = "cluster_6", DOB = "cluster_12", DPL = "cluster_5", DIS = "cluster_13"
  )
)

# =============================================================================
# 1. Probability-weighted regression
# =============================================================================

run_weighted_regression <- function(data, prob_col, biomarker) {
  w <- data[[prob_col]]
  if (all(is.na(w)) || sd(w, na.rm = TRUE) < 1e-10) return(NULL)
  if (!biomarker %in% names(data)) return(NULL)
  tryCatch({
    mod   <- lm(reformulate(c(EXPOSURE, COVARIATES), biomarker),
                data = data, weights = w)
    coefs <- summary(mod)$coefficients
    ci    <- confint(mod)
    if (!EXPOSURE %in% rownames(coefs)) return(NULL)
    data.frame(
      biomarker  = biomarker,
      estimate   = coefs[EXPOSURE, 1],
      std_error  = coefs[EXPOSURE, 2],
      t_value    = coefs[EXPOSURE, 3],
      p_value    = coefs[EXPOSURE, 4],
      conf_low   = ci[EXPOSURE, 1],
      conf_high  = ci[EXPOSURE, 2],
      n_weighted = round(sum(w, na.rm = TRUE), 1)
    )
  }, error = function(e) NULL)
}

cat("=== Cluster-specific TyG effects ===\n")

all_effects <- map_dfr(SEX_GROUPS, function(sg) {
  data    <- complete_results$sex_stratified[[sg]]
  probs   <- complete_results$clustering_results[[sg]]$probs
  mapping <- CLUSTER_MAP[[sg]]

  prob_renamed <- probs |>
    rename_with(~ paste0("prob_", .x), starts_with("cluster_"))
  merged <- inner_join(
    prob_renamed,
    data |> mutate(eid = as.character(eid)),
    by = "eid"
  )

  map_dfr(names(mapping), function(label) {
    prob_col <- paste0("prob_", mapping[[label]])
    if (!prob_col %in% names(merged)) return(NULL)
    map_dfr(BIOMARKERS, ~ run_weighted_regression(merged, prob_col, .x)) |>
      mutate(profile = label)
  }) |>
    mutate(sex = sg)
}) |>
  mutate(is_significant = p_value < 0.05)

cat(sprintf("Completed %d effect estimates.\n", nrow(all_effects)))

# =============================================================================
# 2. Forest plot
# =============================================================================

plot_data <- all_effects |>
  mutate(
    biomarker_display = factor(BIOMARKER_LABELS[biomarker], levels = BIOMARKER_LABELS),
    profile_f         = factor(profile, levels = rev(CLUSTER_ORDER)),
    sex_label         = factor(sex, levels = SEX_GROUPS)
  )

bc_ref <- plot_data |>
  filter(profile == "BC") |>
  select(sex_label, biomarker_display, conf_low, conf_high)

forest_plot <- ggplot(plot_data, aes(x = estimate, y = profile_f)) +
  geom_rect(
    data = bc_ref,
    aes(xmin = conf_low, xmax = conf_high, ymin = -Inf, ymax = Inf),
    fill = "#FFB6C1", alpha = 0.5, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high),
                 linewidth = 1.0, color = "black") +
  geom_point(aes(fill = is_significant), shape = 21,
             size = 3.5, stroke = 0.8, color = "black") +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  facet_grid(sex_label ~ biomarker_display, scales = "free_x", space = "free_y") +
  labs(
    x     = "Difference in biomarker per unit increase in TyG",
    y     = NULL,
    title = "Cluster-specific TyG effects on biomarkers"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(size = 7),
    axis.text.y      = element_text(size = 9, face = "bold"),
    axis.title.x     = element_text(size = 9, margin = margin(t = 6)),
    strip.text.x     = element_text(size = 8, face = "bold"),
    strip.text.y     = element_text(size = 9, face = "bold", angle = 270),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.4),
    plot.title       = element_text(size = 11, hjust = 0.5, face = "bold"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.4),
    panel.spacing    = unit(0.8, "lines")
  )

ggsave(file.path(OUTPUT_DIR, "TyG_cluster_effects.png"),
       forest_plot, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "TyG_cluster_effects.pdf"),
       forest_plot, width = 10, height = 7, bg = "white")

cat("Forest plot saved.\n")

# =============================================================================
# 3. Supplementary Excel table
# =============================================================================

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "f"))
}

BM_ORDER <- toupper(BIOMARKERS)
N_BM     <- length(BM_ORDER)
N_COLS   <- 1 + N_BM * 2

table_wide <- all_effects |>
  mutate(
    ci_str   = sprintf("%.3f (%.3f, %.3f)", estimate, conf_low, conf_high),
    p_str    = format_p(p_value),
    bm_upper = toupper(biomarker)
  ) |>
  select(sex, profile, bm_upper, ci_str, p_str) |>
  pivot_wider(
    names_from  = bm_upper,
    values_from = c(ci_str, p_str),
    names_glue  = "{bm_upper}_{.value}"
  ) |>
  mutate(
    profile = factor(profile, levels = CLUSTER_ORDER),
    sex     = factor(sex, levels = SEX_GROUPS)
  ) |>
  arrange(sex, profile) |>
  mutate(
    profile = paste0(profile, " (", PHENOTYPE_LABELS[as.character(profile)], ")")
  )

wb <- createWorkbook()
addWorksheet(wb, "Cluster TyG Effects")
ws <- "Cluster TyG Effects"

# Cell styles
sty <- list(
  title    = createStyle(fontSize = 11, fontColour = "white", fgFill = "#2F4F6F",
                         halign = "CENTER", textDecoration = "bold", wrapText = TRUE),
  bm_head  = createStyle(fontSize = 9,  fontColour = "white", fgFill = "#5C6B7A",
                         halign = "CENTER", textDecoration = "bold", wrapText = TRUE),
  sub_head = createStyle(fontSize = 8.5, fgFill = "#D9E6F2",
                         halign = "CENTER", textDecoration = "bold"),
  sex_f    = createStyle(fontSize = 10, fontColour = "white", fgFill = "#5B8DB8",
                         halign = "CENTER", textDecoration = "bold"),
  sex_m    = createStyle(fontSize = 10, fontColour = "white", fgFill = "#4A7A6D",
                         halign = "CENTER", textDecoration = "bold"),
  cell     = createStyle(fontSize = 9, halign = "CENTER"),
  cell_sig = createStyle(fontSize = 9, halign = "CENTER",
                         fontColour = "#C0392B", textDecoration = "bold"),
  cell_bc  = createStyle(fontSize = 9, halign = "CENTER", fgFill = "#F2F2F2"),
  cluster  = createStyle(fontSize = 9, halign = "LEFT")
)

# Title row
writeData(wb, ws,
          "Supplementary Table. Cluster-specific TyG\u2013biomarker associations",
          startRow = 1, startCol = 1)
mergeCells(wb, ws, cols = 1:N_COLS, rows = 1)
addStyle(wb, ws, sty$title, rows = 1, cols = 1:N_COLS, gridExpand = TRUE)
setRowHeights(wb, ws, 1, 30)

# Biomarker header row
writeData(wb, ws, "Profile", startRow = 2, startCol = 1)
addStyle(wb, ws, sty$bm_head, rows = 2, cols = 1)
for (i in seq_along(BM_ORDER)) {
  bm  <- tolower(BM_ORDER[i])
  col <- 1 + (i - 1) * 2 + 1
  writeData(wb, ws,
            sprintf("%s\n(%s)", BM_ORDER[i], BIOMARKER_UNITS[bm]),
            startRow = 2, startCol = col)
  mergeCells(wb, ws, cols = col:(col + 1), rows = 2)
  addStyle(wb, ws, sty$bm_head, rows = 2, cols = col:(col + 1), gridExpand = TRUE)
}
setRowHeights(wb, ws, 2, 32)

# Sub-header row
writeData(wb, ws, "Profile", startRow = 3, startCol = 1)
addStyle(wb, ws, sty$sub_head, rows = 3, cols = 1)
for (i in seq_along(BM_ORDER)) {
  col <- 1 + (i - 1) * 2 + 1
  writeData(wb, ws, "\u03b2 (95% CI)", startRow = 3, startCol = col)
  writeData(wb, ws, "p-value",         startRow = 3, startCol = col + 1)
  addStyle(wb, ws, sty$sub_head, rows = 3, cols = col:(col + 1), gridExpand = TRUE)
}
setRowHeights(wb, ws, 3, 22)

# Write a sex block
write_block <- function(df_sex, start_row, sex_label, sex_sty) {
  n <- nrow(df_sex)
  writeData(wb, ws, sex_label, startRow = start_row, startCol = 1)
  mergeCells(wb, ws, cols = 1:N_COLS, rows = start_row)
  addStyle(wb, ws, sex_sty, rows = start_row, cols = 1:N_COLS, gridExpand = TRUE)
  setRowHeights(wb, ws, start_row, 20)

  for (i in seq_len(n)) {
    row_i  <- start_row + i
    is_bc  <- grepl("^BC", df_sex$profile[i])
    writeData(wb, ws, df_sex$profile[i], startRow = row_i, startCol = 1)
    addStyle(wb, ws, sty$cluster, rows = row_i, cols = 1)

    for (j in seq_along(BM_ORDER)) {
      bm     <- BM_ORDER[j]
      col_ci <- 1 + (j - 1) * 2 + 1
      col_p  <- col_ci + 1
      ci_val <- df_sex[[paste0(bm, "_ci_str")]][i]
      p_val  <- df_sex[[paste0(bm, "_p_str")]][i]

      p_num <- all_effects |>
        filter(
          sex      == sub(" \\(.*", "", sex_label),
          profile  == sub(" \\(.*", "", df_sex$profile[i]),
          biomarker == tolower(bm)
        ) |>
        pull(p_value)
      is_sig <- length(p_num) > 0 && !is.na(p_num[1]) && p_num[1] < 0.05

      cell_sty <- if (is_bc) sty$cell_bc else if (is_sig) sty$cell_sig else sty$cell
      writeData(wb, ws, ci_val, startRow = row_i, startCol = col_ci)
      writeData(wb, ws, p_val,  startRow = row_i, startCol = col_p)
      addStyle(wb, ws, cell_sty, rows = row_i, cols = col_ci:col_p, gridExpand = TRUE)
    }
    setRowHeights(wb, ws, row_i, 18)
  }
  start_row + n
}

last_f <- write_block(filter(table_wide, sex == "Female"), 4, "Female", sty$sex_f)
last_m <- write_block(filter(table_wide, sex == "Male"), last_f + 2, "Male", sty$sex_m)

# Footnote
note_row <- last_m + 2
writeData(wb, ws,
  paste(
    "Note: Values are regression coefficients \u03b2 (95% CI) from",
    "probability-weighted linear regression within each TyG-discordant profile.",
    "Bold red: p < 0.05 (uncorrected). Grey shading: BC (reference).",
    "Covariates: age, smoking."
  ),
  startRow = note_row, startCol = 1
)
mergeCells(wb, ws, cols = 1:N_COLS, rows = note_row)
addStyle(wb, ws,
         createStyle(fontSize = 8, fontColour = "#555555", wrapText = TRUE),
         rows = note_row, cols = 1:N_COLS, gridExpand = TRUE)
setRowHeights(wb, ws, note_row, 60)

# Column widths and freeze pane
setColWidths(wb, ws, cols = 1, widths = 36)
for (i in seq_along(BM_ORDER)) {
  col <- 1 + (i - 1) * 2 + 1
  setColWidths(wb, ws, cols = col,     widths = 22)
  setColWidths(wb, ws, cols = col + 1, widths = 9)
}
freezePane(wb, ws, firstActiveRow = 5, firstActiveCol = 2)

out_xlsx <- file.path(OUTPUT_DIR, "Supplementary_Table_Cluster_TyG_Effects.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

write.csv(all_effects,
          file.path(OUTPUT_DIR, "TyG_cluster_effects.csv"),
          row.names = FALSE)

cat("Done. Saved: TyG_cluster_effects.png/.pdf, .csv, and .xlsx\n")
