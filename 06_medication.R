# =============================================================================
# 06_medication.R
# TyG-Biomarker Discordance Profiling
# Step 6: Cluster × medication use associations
#         ALR-logistic regression + LRT interaction test + forest plot
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
#           imp_data1.csv
# Output  : TyG_medication_forest.png / .pdf
#           TyG_medication_results.csv
#           TyG_medication_LRT.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(scales)

OUTPUT_DIR <- "."
DATA_PATH  <- "imp_data1.csv"

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))
imp_data1 <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# =============================================================================
# Constants
# =============================================================================

MEDICATIONS <- c("lipid_lowering", "antidiabetic_total", "antihypertensive")

MED_LABELS <- c(
  lipid_lowering     = "Lipid-Lowering Therapy",
  antidiabetic_total = "Antidiabetic Therapy",
  antihypertensive   = "Antihypertensive Therapy"
)

MED_ORDER <- unname(MED_LABELS)

MED_SHORT <- c(
  "Lipid-Lowering Therapy"   = "LipidLower",
  "Antidiabetic Therapy"     = "AntiDM",
  "Antihypertensive Therapy" = "AntiHT"
)

COVARIATES_ADJ <- c("age", "smoking_status", "drinking.status", "TDI", "Education_Level")
CLUSTER_ORDER  <- c("DHT", "DRI", "DHL", "DOB", "DPL", "DIS")
SEX_GROUPS     <- c("Female", "Male")
MIN_CASES      <- 5

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
# 1. Data integration and ALR transformation
# =============================================================================

prepare_sex_data <- function(sex_name) {
  profile <- CLUSTER_MAP[[sex_name]]
  probs   <- complete_results$clustering_results[[sex_name]]$probs |>
    select(eid, all_of(unname(unlist(profile)))) |>
    rename(all_of(setNames(unname(unlist(profile)), names(unlist(profile)))))

  sex_values <- if (sex_name == "Female") c("F", "Female", "female", "0") else
    c("M", "Male", "male", "1")

  pheno <- imp_data1 |>
    filter(sex %in% sex_values) |>
    select(eid, TyG, all_of(c(MEDICATIONS, COVARIATES_ADJ)))

  inner_join(probs, pheno, by = "eid")
}

do_alr <- function(data, cl_order = CLUSTER_ORDER) {
  all_cls  <- c("BC", cl_order)
  prob_mat <- pmax(pmin(as.matrix(data[, all_cls]), 1 - 1e-15), 1e-15)
  alr_mat  <- log(prob_mat[, cl_order, drop = FALSE] / prob_mat[, "BC"])
  colnames(alr_mat) <- paste0("alr_", cl_order)
  bind_cols(select(data, -all_of(all_cls)), as.data.frame(alr_mat))
}

cat("=== Data integration ===\n")
data_female_alr <- do_alr(prepare_sex_data("Female"))
data_male_alr   <- do_alr(prepare_sex_data("Male"))

# =============================================================================
# 2. Logistic regression + LRT interaction test
# =============================================================================

run_logistic_med <- function(data_alr, medication) {
  alr_terms <- paste0("alr_", CLUSTER_ORDER)
  n_cases   <- sum(data_alr[[medication]], na.rm = TRUE)

  if (n_cases <= MIN_CASES) {
    message(sprintf("  Skipping %s (n = %d)", medication, n_cases))
    return(list(results = NULL, lrt = NULL))
  }

  extract_glm <- function(mod, label) {
    if (is.null(mod)) return(NULL)
    coefs <- summary(mod)$coefficients
    rows  <- grep("^alr_", rownames(coefs))
    if (!length(rows)) return(NULL)
    est <- coefs[rows, "Estimate"]
    se  <- coefs[rows, "Std. Error"]
    data.frame(
      subtype          = gsub("^alr_", "", rownames(coefs)[rows]),
      log_or = est, se = se, pval = coefs[rows, "Pr(>|z|)"],
      or = exp(est), or_lower = exp(est - qnorm(0.975) * se),
      or_upper = exp(est + qnorm(0.975) * se),
      model = label, medication = medication,
      medication_label = MED_LABELS[medication],
      n_users = n_cases
    )
  }

  mod_crude <- tryCatch(
    glm(reformulate(alr_terms, medication), data_alr, family = binomial()),
    error = function(e) NULL
  )
  mod_adj <- tryCatch(
    glm(reformulate(c(alr_terms, COVARIATES_ADJ), medication),
        data_alr, family = binomial()),
    error = function(e) NULL
  )

  # LRT
  lrt_row <- tibble(
    medication = medication, n_lrt = NA_real_,
    lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_
  )

  if ("TyG" %in% names(data_alr)) {
    data_lrt <- data_alr |>
      mutate(TyG_c = TyG - mean(TyG, na.rm = TRUE)) |>
      select(all_of(c(medication, alr_terms, COVARIATES_ADJ, "TyG_c"))) |>
      drop_na()

    mod_adj_lrt <- tryCatch(
      glm(reformulate(c(alr_terms, COVARIATES_ADJ), medication),
          data_lrt, family = binomial()),
      error = function(e) NULL
    )
    mod_int <- tryCatch(
      glm(reformulate(c(alr_terms, "TyG_c", COVARIATES_ADJ,
                        paste0(alr_terms, ":TyG_c")), medication),
          data_lrt, family = binomial()),
      error = function(e) NULL
    )

    lrt_row$n_lrt <- nrow(data_lrt)
    if (!is.null(mod_adj_lrt) && !is.null(mod_int)) {
      lrt <- tryCatch(anova(mod_adj_lrt, mod_int, test = "Chisq"), error = function(e) NULL)
      if (!is.null(lrt) && nrow(lrt) >= 2) {
        lrt_row <- tibble(
          medication = medication,     n_lrt     = nrow(data_lrt),
          lrt_chisq  = lrt[["Deviance"]][2],
          lrt_df     = lrt[["Df"]][2],
          lrt_p      = lrt[["Pr(>Chi))"]][2]
        )
      }
    }
  }

  list(
    results = bind_rows(extract_glm(mod_crude, "Unadjusted"),
                        extract_glm(mod_adj,   "Adjusted")),
    lrt = lrt_row
  )
}

# =============================================================================
# 3. Run regressions
# =============================================================================

cat("\n=== Logistic regression ===\n")

run_all <- function(data_alr, sex_name) {
  out <- map(MEDICATIONS, ~ run_logistic_med(data_alr, .x))
  list(
    results = map_dfr(out, "results") |> mutate(sex = sex_name),
    lrt     = map_dfr(out, "lrt")     |> mutate(sex = sex_name)
  )
}

res_f <- run_all(data_female_alr, "Female")
res_m <- run_all(data_male_alr,   "Male")

med_col_order <- unname(MED_SHORT[MED_ORDER])

results_all <- bind_rows(res_f$results, res_m$results) |>
  mutate(
    subtype          = factor(subtype, levels = CLUSTER_ORDER),
    medication_label = factor(medication_label, levels = MED_ORDER),
    model            = factor(model, levels = c("Unadjusted", "Adjusted")),
    sex              = factor(sex, levels = SEX_GROUPS),
    med_short        = MED_SHORT[as.character(medication_label)],
    med_short        = factor(med_short, levels = unname(MED_SHORT[MED_ORDER]))
  ) |>
  group_by(sex, model, medication) |>
  mutate(p_fdr = p.adjust(pval, method = "BH"), sig_fdr = p_fdr < 0.05) |>
  ungroup()

lrt_results <- bind_rows(res_f$lrt, res_m$lrt) |>
  mutate(medication_label = factor(MED_LABELS[medication], levels = MED_ORDER)) |>
  group_by(sex) |>
  mutate(lrt_p_fdr = p.adjust(lrt_p, method = "BH")) |>
  ungroup()

cat("\n=== LRT interaction test results ===\n")
print(select(lrt_results, sex, medication_label, lrt_chisq, lrt_df, lrt_p, lrt_p_fdr))

# =============================================================================
# 4. Forest plot
# =============================================================================

MODEL_COLORS <- c("Unadjusted" = "black", "Adjusted" = "#E41A1C")

x_axis_settings <- list(
  LipidLower = list(limits = c(0.75, 1.30), breaks = c(0.80, 1.00, 1.20)),
  AntiDM     = list(limits = c(0.75, 1.30), breaks = c(0.80, 1.00, 1.20)),
  AntiHT     = list(limits = c(0.85, 1.15), breaks = c(0.90, 1.00, 1.10))
)

plot_data <- results_all |>
  rowwise() |>
  mutate(
    lim_lo        = x_axis_settings[[as.character(med_short)]]$limits[1],
    lim_hi        = x_axis_settings[[as.character(med_short)]]$limits[2],
    or_lower_plot = pmax(or_lower, lim_lo),
    or_upper_plot = pmin(or_upper, lim_hi)
  ) |>
  ungroup() |>
  mutate(subtype = factor(subtype, levels = rev(CLUSTER_ORDER)))

x_scales_list <- map(med_col_order, function(m) {
  s <- x_axis_settings[[m]]
  scale_x_log10(limits = s$limits, breaks = s$breaks,
                labels = number_format(accuracy = 0.01))
}) |> set_names(med_col_order)

forest_plot <- ggplot(
  plot_data,
  aes(x = or, y = subtype, color = model, group = model, shape = sig_fdr)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = or_lower_plot, xmax = or_upper_plot),
                 linewidth = 0.8, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = MODEL_COLORS,
                     name   = NULL,
                     breaks = c("Adjusted", "Unadjusted")) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "FDR q < 0.05", "FALSE" = "FDR q \u2265 0.05"),
    name   = "Significance"
  ) +
  facet_grid2(sex ~ med_short, scales = "free_x") +
  facetted_pos_scales(x = x_scales_list) +
  labs(x = "OR relative to BC", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    strip.text.x     = element_text(size = 8.5, face = "bold"),
    strip.text.y     = element_text(size = 9,   face = "bold", angle = 270),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.4),
    axis.text.x      = element_text(size = 7),
    axis.title.x     = element_text(size = 8.5, margin = margin(t = 8)),
    axis.text.y      = element_text(size = 8, face = "bold"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.spacing    = unit(0.8, "lines"),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_text(size = 8.5),
    plot.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16, size = 3)),
    shape = guide_legend(title = "Significance")
  )

ggsave(file.path(OUTPUT_DIR, "TyG_medication_forest.png"),
       forest_plot, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "TyG_medication_forest.pdf"),
       forest_plot, width = 10, height = 7, bg = "white")

write.csv(results_all, file.path(OUTPUT_DIR, "TyG_medication_results.csv"), row.names = FALSE)
write.csv(lrt_results, file.path(OUTPUT_DIR, "TyG_medication_LRT.csv"),     row.names = FALSE)

cat("Done. All files saved to:", OUTPUT_DIR, "\n")
