# =============================================================================
# 05_prevalent_disease.R
# TyG-Biomarker Discordance Profiling
# Step 5: Cluster × prevalent disease associations
#         ALR-logistic regression + LRT interaction test + forest plot
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
#           imp_data1.csv
# Output  : TyG_prevalent_disease_forest.png / .pdf
#           TyG_prevalent_disease_results.csv
#           TyG_prevalent_disease_LRT.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)

OUTPUT_DIR <- "."
DATA_PATH  <- "imp_data1.csv"

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))
imp_data1 <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# =============================================================================
# Constants
# =============================================================================

DISEASES <- c(
  "D1M_baseline", "D2M_baseline", "HT_baseline",
  "RA_baseline",  "stroke_baseline", "CHD_baseline"
)

DISEASE_LABELS <- c(
  D1M_baseline    = "Type 1 Diabetes",
  D2M_baseline    = "Type 2 Diabetes",
  HT_baseline     = "Hypertension",
  RA_baseline     = "Rheumatoid Arthritis",
  stroke_baseline = "Stroke",
  CHD_baseline    = "Coronary Heart Disease"
)

MEDICATIONS    <- c("lipid_lowering", "antidiabetic_total", "antihypertensive")
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
# 1. Data integration
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
    select(eid, TyG, all_of(c(DISEASES, MEDICATIONS, COVARIATES_ADJ)))

  joined <- inner_join(probs, pheno, by = "eid")
  cat(sprintf("  %s: n = %d\n", sex_name, nrow(joined)))
  joined
}

cat("=== Data integration ===\n")
data_female <- prepare_sex_data("Female")
data_male   <- prepare_sex_data("Male")

# =============================================================================
# 2. Additive log-ratio (ALR) transformation
# =============================================================================

do_alr <- function(data, cl_order = CLUSTER_ORDER) {
  all_cls  <- c("BC", cl_order)
  prob_mat <- as.matrix(data[, all_cls])
  prob_mat <- pmax(pmin(prob_mat, 1 - 1e-15), 1e-15)
  alr_mat  <- log(prob_mat[, cl_order, drop = FALSE] / prob_mat[, "BC"])
  colnames(alr_mat) <- paste0("alr_", cl_order)
  bind_cols(select(data, -all_of(all_cls)), as.data.frame(alr_mat))
}

data_female_alr <- do_alr(data_female)
data_male_alr   <- do_alr(data_male)

# =============================================================================
# 3. Logistic regression + LRT interaction test
# =============================================================================

extract_glm_coefs <- function(mod, model_label, disease, n_cases) {
  if (is.null(mod)) return(NULL)
  coefs <- summary(mod)$coefficients
  rows  <- grep("^alr_", rownames(coefs))
  if (!length(rows)) return(NULL)
  est <- coefs[rows, "Estimate"]
  se  <- coefs[rows, "Std. Error"]
  data.frame(
    subtype  = gsub("^alr_", "", rownames(coefs)[rows]),
    log_or   = est, se = se,
    pval     = coefs[rows, "Pr(>|z|)"],
    or       = exp(est),
    or_lower = exp(est - qnorm(0.975) * se),
    or_upper = exp(est + qnorm(0.975) * se),
    model    = model_label,
    disease  = disease,
    n_cases  = n_cases
  )
}

run_logistic <- function(data_alr, disease) {
  alr_terms <- paste0("alr_", CLUSTER_ORDER)
  n_cases   <- sum(data_alr[[disease]], na.rm = TRUE)

  if (n_cases <= MIN_CASES) {
    message(sprintf("  Skipping %s (n cases = %d)", disease, n_cases))
    return(list(results = NULL, lrt = NULL))
  }

  mod_crude <- tryCatch(
    glm(reformulate(alr_terms, disease), data_alr, family = binomial()),
    error = function(e) NULL
  )
  mod_adj <- tryCatch(
    glm(reformulate(c(alr_terms, COVARIATES_ADJ, MEDICATIONS), disease),
        data_alr, family = binomial()),
    error = function(e) NULL
  )

  # LRT: test cluster × TyG interaction
  data_lrt <- data_alr |>
    mutate(TyG_c = TyG - mean(TyG, na.rm = TRUE)) |>
    select(all_of(c(disease, alr_terms, COVARIATES_ADJ, MEDICATIONS, "TyG_c"))) |>
    drop_na()

  lrt_row <- tibble(
    disease   = disease, n_lrt = nrow(data_lrt),
    lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_
  )

  mod_adj_lrt <- tryCatch(
    glm(reformulate(c(alr_terms, COVARIATES_ADJ, MEDICATIONS, "TyG_c"), disease),
        data_lrt, family = binomial()),
    error = function(e) NULL
  )
  mod_int <- tryCatch(
    glm(reformulate(c(alr_terms, "TyG_c", COVARIATES_ADJ, MEDICATIONS,
                      paste0(alr_terms, ":TyG_c")), disease),
        data_lrt, family = binomial()),
    error = function(e) NULL
  )

  if (!is.null(mod_adj_lrt) && !is.null(mod_int)) {
    lrt <- tryCatch(anova(mod_adj_lrt, mod_int, test = "Chisq"), error = function(e) NULL)
    if (!is.null(lrt) && nrow(lrt) >= 2) {
      lrt_row <- tibble(
        disease   = disease,     n_lrt     = nrow(data_lrt),
        lrt_chisq = lrt[["Deviance"]][2],
        lrt_df    = lrt[["Df"]][2],
        lrt_p     = lrt[["Pr(>Chi))"]][2]
      )
    }
  }

  list(
    results = bind_rows(
      extract_glm_coefs(mod_crude, "Unadjusted", disease, n_cases),
      extract_glm_coefs(mod_adj,   "Adjusted",   disease, n_cases)
    ),
    lrt = lrt_row
  )
}

# =============================================================================
# 4. Run regressions
# =============================================================================

cat("\n=== Logistic regression ===\n")

run_all <- function(data_alr, sex_name) {
  out <- map(DISEASES, ~ run_logistic(data_alr, .x))
  list(
    results = map_dfr(out, "results") |> mutate(sex = sex_name),
    lrt     = map_dfr(out, "lrt")     |> mutate(sex = sex_name)
  )
}

res_f <- run_all(data_female_alr, "Female")
res_m <- run_all(data_male_alr,   "Male")

results_all <- bind_rows(res_f$results, res_m$results) |>
  mutate(
    subtype       = factor(subtype, levels = CLUSTER_ORDER),
    disease_label = factor(DISEASE_LABELS[disease], levels = unname(DISEASE_LABELS)),
    model         = factor(model, levels = c("Unadjusted", "Adjusted"))
  ) |>
  group_by(sex, model, disease) |>
  mutate(p_fdr = p.adjust(pval, method = "BH"), sig_fdr = p_fdr < 0.05) |>
  ungroup()

lrt_results <- bind_rows(res_f$lrt, res_m$lrt) |>
  mutate(disease_label = factor(DISEASE_LABELS[disease], levels = unname(DISEASE_LABELS))) |>
  group_by(sex) |>
  mutate(lrt_p_fdr = p.adjust(lrt_p, method = "BH")) |>
  ungroup()

cat("\n=== LRT interaction test results ===\n")
print(select(lrt_results, sex, disease_label, lrt_chisq, lrt_df, lrt_p, lrt_p_fdr))

# =============================================================================
# 5. Forest plot
# =============================================================================

MODEL_COLORS <- c("Unadjusted" = "grey50", "Adjusted" = "#E41A1C")

x_lim <- results_all |>
  filter(is.finite(or_lower), is.finite(or_upper)) |>
  summarise(lo = quantile(or_lower, 0.01, na.rm = TRUE),
            hi = quantile(or_upper, 0.99, na.rm = TRUE))
X_LIM <- c(max(0.05, x_lim$lo * 0.8), min(50, x_lim$hi * 1.2))

make_panel <- function(sex_name, disease_name) {
  viz <- results_all |>
    filter(sex == sex_name, disease_label == disease_name) |>
    mutate(
      subtype      = factor(subtype, levels = rev(CLUSTER_ORDER)),
      or_lower_plt = pmax(or_lower, X_LIM[1]),
      or_upper_plt = pmin(or_upper, X_LIM[2])
    )
  if (!nrow(viz)) return(ggplot() + theme_void())

  n_lab <- viz |>
    distinct(disease, n_cases) |>
    mutate(label = paste0("n cases = ", n_cases))

  ggplot(viz, aes(y = subtype, x = or, color = model, shape = sig_fdr)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = or_lower_plt, xmax = or_upper_plt),
                   height = 0.15, linewidth = 0.45,
                   position = position_dodge(width = 0.65)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.65)) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
    scale_color_manual(values = MODEL_COLORS, name = "Model") +
    scale_x_log10(limits = X_LIM) +
    geom_text(data = n_lab,
              aes(x = X_LIM[2], y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 1.5,
              size = 2.2, color = "grey40") +
    labs(x = "OR (vs BC)", y = NULL, title = disease_name) +
    theme_bw(base_size = 9) +
    theme(
      plot.title      = element_text(size = 9, hjust = 0.5),
      axis.text.x     = element_text(size = 7),
      axis.text.y     = element_text(size = 8),
      panel.grid      = element_blank(),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
      legend.position = "none",
      plot.margin     = margin(2, 3, 2, 3, "mm")
    )
}

build_row <- function(sex_name) {
  panels <- map(unname(DISEASE_LABELS), ~ make_panel(sex_name, .x))
  sex_strip <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = sex_name,
             angle = 270, size = 3.5, fontface = "bold") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey90", color = NA))
  do.call(plot_grid, c(panels, list(sex_strip),
                       list(nrow = 1,
                            rel_widths = c(rep(1, length(DISEASES)), 0.1),
                            align = "h", axis = "tb")))
}

legend_plot <- ggplot(
  data.frame(
    model = factor(rep(c("Unadjusted", "Adjusted"), 2),
                   levels = c("Unadjusted", "Adjusted")),
    x = 1, y = 1:4
  ),
  aes(x, y, color = model)
) +
  geom_point(size = 2.5) +
  scale_color_manual(values = MODEL_COLORS, name = "Model") +
  theme_void() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 8.5))

final_plot <- plot_grid(
  get_legend(legend_plot),
  build_row("Female"),
  build_row("Male"),
  ncol = 1, rel_heights = c(0.07, 1, 1)
)

ggsave(file.path(OUTPUT_DIR, "TyG_prevalent_disease_forest.png"),
       final_plot, width = 18, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "TyG_prevalent_disease_forest.pdf"),
       final_plot, width = 18, height = 8)

# =============================================================================
# 6. Save results
# =============================================================================

write.csv(results_all,
          file.path(OUTPUT_DIR, "TyG_prevalent_disease_results.csv"),
          row.names = FALSE)
write.csv(lrt_results,
          file.path(OUTPUT_DIR, "TyG_prevalent_disease_LRT.csv"),
          row.names = FALSE)

cat("Done. All files saved to:", OUTPUT_DIR, "\n")
