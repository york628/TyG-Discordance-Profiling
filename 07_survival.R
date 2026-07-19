# =============================================================================
# 07_survival.R
# TyG-Biomarker Discordance Profiling
# Step 7: Cluster × MACE and incident T2D survival analysis
#         ALR-Cox models + forest plots
#
# Requires: TyG_clustering_results.RData (from 01_clustering.R)
#           imp_data1.csv
# Output  : TyG_survival_forest.png / .pdf
#           TyG_survival_HR.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(ggh4x)
library(scales)

OUTPUT_DIR   <- "."
DATA_PATH    <- "imp_data1.csv"
STUDY_CUTOFF <- as.Date("2024-11-30")
FUT_YEAR     <- 10   # follow-up horizon (years)

load(file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))
imp_data1 <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# =============================================================================
# Constants
# =============================================================================

CLUSTER_ORDER <- c("DHT", "DRI", "DHL", "DOB", "DPL", "DIS")
SEX_GROUPS    <- c("Female", "Male")

COVARS_MACE <- c(
  "age", "smoking_current", "drinking.status", "TDI", "Education_Level",
  "lipid_lowering", "antihypertensive", "antidiabetic_total"
)
COVARS_T2D <- c(
  "age", "smoking_current", "drinking.status", "TDI", "Education_Level",
  "antihypertensive"
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
# 1. Outcome construction
# =============================================================================

cat("=== Constructing outcomes ===\n")

imp_data1 <- imp_data1 |>
  mutate(
    smoking_current = if_else(
      smoking_status %in% c("Only occasionally", "Yes, on most or all days"),
      1L, 0L, missing = 0L
    ),
    date_baseline = as.Date(Date.of.attending.assessment.centre0)
  )

build_mace <- function(df) {
  df |>
    mutate(
      date_death  = as.Date(death),
      date_ltfu   = as.Date(Date.lost.to.follow.up),
      date_stroke = as.Date(stroke_time),
      date_MI     = as.Date(MI_time),
      mace_occurred = (fatal_mace == 1L | stroke == 1L | MI == 1L),
      date_mace_event = pmin(
        if_else(fatal_mace == 1L, date_death,  as.Date(NA)),
        if_else(stroke     == 1L, date_stroke, as.Date(NA)),
        if_else(MI         == 1L, date_MI,     as.Date(NA)),
        na.rm = TRUE
      ),
      date_end = pmin(
        if_else(mace_occurred, date_mace_event, as.Date(NA)),
        if_else(!mace_occurred & !is.na(date_death), date_death, as.Date(NA)),
        date_ltfu, STUDY_CUTOFF, na.rm = TRUE
      ),
      mace_timeyrs = as.numeric(date_end - date_baseline) / 365.25,
      mace_event   = as.integer(
        mace_occurred & !is.na(date_mace_event) & date_mace_event == date_end
      )
    ) |>
    filter(!is.na(mace_timeyrs), mace_timeyrs > 0)
}

build_t2d <- function(df) {
  df |>
    mutate(
      date_ltfu   = as.Date(Date.lost.to.follow.up),
      date_t2d    = as.Date(D2M_time),
      date_end    = pmin(if_else(D2M == 1L, date_t2d, as.Date(NA)),
                         date_ltfu, STUDY_CUTOFF, na.rm = TRUE),
      t2d_timeyrs = as.numeric(date_end - date_baseline) / 365.25,
      t2d_event   = as.integer(D2M == 1L & !is.na(date_t2d) & date_t2d == date_end)
    ) |>
    filter(D1M_baseline == 0, D2M_baseline == 0,
           !is.na(t2d_timeyrs), t2d_timeyrs > 0)
}

imp_data1 <- build_mace(imp_data1) |> build_t2d()

cat(sprintf("MACE events: %d\nT2D  events: %d\n",
            sum(imp_data1$mace_event, na.rm = TRUE),
            sum(imp_data1$t2d_event,  na.rm = TRUE)))

# =============================================================================
# 2. Data integration and ALR transformation
# =============================================================================

prepare_sex_data <- function(sex_name) {
  profile <- CLUSTER_MAP[[sex_name]]
  probs   <- complete_results$clustering_results[[sex_name]]$probs |>
    select(eid, all_of(unname(unlist(profile)))) |>
    rename(all_of(setNames(unname(unlist(profile)), names(unlist(profile)))))

  sex_values <- if (sex_name == "Female") c("F", "Female", "female", "0") else
    c("M", "Male", "male", "1")

  inner_join(probs,
             filter(imp_data1, as.character(sex) %in% sex_values),
             by = "eid")
}

do_alr <- function(data, cl_order = CLUSTER_ORDER) {
  all_cls  <- c("BC", cl_order)
  prob_mat <- pmax(pmin(as.matrix(data[, all_cls]), 1 - 1e-15), 1e-15)
  alr_mat  <- log(prob_mat[, cl_order, drop = FALSE] / prob_mat[, "BC"])
  colnames(alr_mat) <- paste0("alr_", cl_order)
  bind_cols(select(data, -all_of(all_cls)), as.data.frame(alr_mat))
}

data_female <- prepare_sex_data("Female")
data_male   <- prepare_sex_data("Male")

# =============================================================================
# 3. Cox model helpers
# =============================================================================

fit_cox <- function(formula_str, data) {
  tryCatch(
    coxph(as.formula(formula_str), data = data),
    error = function(e) { message("  [Cox error] ", e$message); NULL }
  )
}

extract_alr_hr <- function(mod, type_label) {
  if (is.null(mod)) return(NULL)
  coefs <- coef(mod)
  vc    <- vcov(mod)
  nms   <- names(coefs)[startsWith(names(coefs), "alr_")]
  map_dfr(nms, function(nm) {
    log_hr <- coefs[nm]; se <- sqrt(vc[nm, nm])
    tibble(
      type    = type_label,
      Cluster = gsub("^alr_", "", nm),
      HR      = exp(log_hr),
      HR_lo   = exp(log_hr - 1.96 * se),
      HR_hi   = exp(log_hr + 1.96 * se),
      p.value = 2 * pnorm(-abs(log_hr / se))
    )
  })
}

calc_cox_hr <- function(data_full, outcome_time, outcome_event, covars) {
  alr_cols <- paste0("alr_", CLUSTER_ORDER)
  dat_alr  <- do_alr(data_full) |>
    select(all_of(intersect(
      c(outcome_time, outcome_event, covars, alr_cols), names(.)
    ))) |>
    mutate(across(everything(), ~ { .x[is.infinite(.x)] <- NaN; .x })) |>
    drop_na() |>
    mutate(
      !!outcome_event := .data[[outcome_event]] * (.data[[outcome_time]] <= FUT_YEAR),
      !!outcome_time  := pmin(.data[[outcome_time]], FUT_YEAR)
    )

  N      <- nrow(dat_alr)
  Ncases <- sum(dat_alr[[outcome_event]])
  cat(sprintf("    N = %d, events = %d\n", N, Ncases))
  if (N < 50 || Ncases < 5) { warning("Insufficient sample size."); return(NULL) }

  surv_str  <- sprintf("Surv(%s, %s)", outcome_time, outcome_event)
  adj_cols  <- intersect(covars, names(dat_alr))
  alr_str   <- paste(alr_cols, collapse = " + ")
  adj_str   <- paste(c(alr_cols, adj_cols), collapse = " + ")

  mod_crude <- fit_cox(paste(surv_str, "~", alr_str),   dat_alr)
  mod_adj   <- fit_cox(paste(surv_str, "~", adj_str),   dat_alr)

  bind_rows(
    extract_alr_hr(mod_crude, "Crude HR"),
    extract_alr_hr(mod_adj,   "Adjusted HR")
  ) |>
    mutate(N = N, Ncases = Ncases)
}

# =============================================================================
# 4. Run Cox models
# =============================================================================

cat("\n=== Cox survival analysis ===\n")

cat("  MACE - Female\n")
hr_mace_f <- calc_cox_hr(
  filter(data_female, CHD_baseline == 0, stroke_baseline == 0),
  "mace_timeyrs", "mace_event", COVARS_MACE
) |> mutate(sex = "Female", outcome = "MACE")

cat("  MACE - Male\n")
hr_mace_m <- calc_cox_hr(
  filter(data_male, CHD_baseline == 0, stroke_baseline == 0),
  "mace_timeyrs", "mace_event", COVARS_MACE
) |> mutate(sex = "Male", outcome = "MACE")

cat("  T2D - Female\n")
hr_t2d_f <- calc_cox_hr(data_female, "t2d_timeyrs", "t2d_event", COVARS_T2D) |>
  mutate(sex = "Female", outcome = "T2D")

cat("  T2D - Male\n")
hr_t2d_m <- calc_cox_hr(data_male, "t2d_timeyrs", "t2d_event", COVARS_T2D) |>
  mutate(sex = "Male", outcome = "T2D")

hr_results <- bind_rows(hr_mace_f, hr_mace_m, hr_t2d_f, hr_t2d_m) |>
  mutate(
    Cluster = factor(Cluster, levels = CLUSTER_ORDER),
    type    = factor(type,    levels = c("Crude HR", "Adjusted HR")),
    sex     = factor(sex,     levels = SEX_GROUPS),
    outcome = factor(outcome, levels = c("MACE", "T2D"))
  ) |>
  group_by(sex, type, outcome) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH"), sig_fdr = p_fdr < 0.05) |>
  ungroup()

cat("\n=== Adjusted HR summary ===\n")
hr_results |>
  filter(type == "Adjusted HR") |>
  mutate(HR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", HR, HR_lo, HR_hi)) |>
  select(outcome, sex, Cluster, HR_fmt, p_fdr, sig_fdr) |>
  print(n = Inf)

# =============================================================================
# 5. Forest plot
# =============================================================================

x_settings <- list(
  MACE = list(limits = c(0.93, 1.06), breaks = c(0.95, 1.00, 1.05)),
  T2D  = list(limits = c(0.89, 1.13), breaks = c(0.90, 0.95, 1.00, 1.05, 1.10))
)

plot_data <- hr_results |>
  mutate(outcome_char = as.character(outcome)) |>
  rowwise() |>
  mutate(
    lim_lo   = x_settings[[outcome_char]]$limits[1],
    lim_hi   = x_settings[[outcome_char]]$limits[2],
    HR_lo_pl = pmax(HR_lo, lim_lo),
    HR_hi_pl = pmin(HR_hi, lim_hi)
  ) |>
  ungroup() |>
  mutate(Cluster = factor(Cluster, levels = rev(CLUSTER_ORDER)))

x_scales_list <- map(c("MACE", "T2D"), function(oc) {
  s <- x_settings[[oc]]
  scale_x_log10(limits = s$limits, breaks = s$breaks,
                labels = number_format(accuracy = 0.01))
}) |> set_names(c("MACE", "T2D"))

MODEL_COLORS <- c("Crude HR" = "black", "Adjusted HR" = "#E41A1C")

forest_plot <- ggplot(
  plot_data,
  aes(x = HR, y = Cluster, color = type, group = type, shape = sig_fdr)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = HR_lo_pl, xmax = HR_hi_pl),
                 linewidth = 0.8, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = MODEL_COLORS,
    labels = c("Crude HR" = "Unadjusted", "Adjusted HR" = "Adjusted"),
    breaks = c("Adjusted HR", "Crude HR"),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "FDR q < 0.05", "FALSE" = "FDR q \u2265 0.05"),
    name   = "Significance"
  ) +
  facet_grid2(sex ~ outcome, scales = "free_x") +
  facetted_pos_scales(x = x_scales_list) +
  labs(x = "HR relative to BC", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    strip.text.x     = element_text(size = 8.5, face = "bold"),
    strip.text.y     = element_text(size = 9, face = "bold", angle = 270),
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
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(OUTPUT_DIR, "TyG_survival_forest.png"),
       forest_plot, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "TyG_survival_forest.pdf"),
       forest_plot, width = 10, height = 7, bg = "white")

write.csv(
  hr_results |> mutate(across(c(HR, HR_lo, HR_hi, p.value, p_fdr), ~ round(.x, 4))),
  file.path(OUTPUT_DIR, "TyG_survival_HR.csv"),
  row.names = FALSE
)

cat("Done. All files saved to:", OUTPUT_DIR, "\n")
