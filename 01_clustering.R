# =============================================================================
# 01_clustering.R
# TyG-Biomarker Discordance Profiling
# Step 1: Data preparation · UMAP clustering · GMM soft assignment
#
# Output: TyG_clustering_results.RData
# =============================================================================

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(readr)
library(uwot)
library(igraph)
library(mvtnorm)

# =============================================================================
# 0. Paths and global constants
# =============================================================================

DATA_PATH   <- "imp_data1.csv"
OUTPUT_DIR  <- "."
RANDOM_SEED <- 42

ID_COL      <- "eid"
SEX_VAR     <- "sex"
SEX_GROUPS  <- c("Female", "Male")
EXPOSURE    <- "tyg"
COVARIATES  <- c("age", "smoking", "drinking", "TDI", "education")
BIOMARKERS  <- c("bmi", "sbp", "dbp", "egfr", "crp", "hdl", "ldl", "chol")
ALL_VARS    <- c(ID_COL, SEX_VAR, EXPOSURE, COVARIATES, BIOMARKERS)

# =============================================================================
# 1. Helper functions
# =============================================================================

# CKD-EPI 2021 eGFR (race-free; serum creatinine in mg/dL)
calc_egfr <- function(scr, age, sex) {
  kappa      <- ifelse(sex == "Female", 0.7,   0.9)
  alpha      <- ifelse(sex == "Female", -0.241, -0.302)
  sex_factor <- ifelse(sex == "Female", 1.012,  1.000)
  ratio      <- scr / kappa
  142 *
    ifelse(ratio < 1, ratio^alpha, 1) *
    ifelse(ratio > 1, ratio^(-1.200), 1) *
    (0.9938^age) * sex_factor
}

# Winsorise at ±5 SD (values outside the window become NaN → treated as missing)
remove_outliers <- function(x, sd_units = 5) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x,   na.rm = TRUE)
  ifelse(x > m - sd_units * s & x < m + sd_units * s, x, NaN)
}

# Extract coefficients, SE, t, p, 95 % CI from a fitted lm
extract_coefs <- function(mod) {
  sm <- summary(mod)$coefficients[, 1:4, drop = FALSE]
  ci <- confint(mod)[, 1:2, drop = FALSE]
  out <- round(cbind(sm, ci), 5)
  colnames(out) <- c("Estimate", "SE", "t_value", "p_value", "lowerCI", "upperCI")
  data.frame(term = rownames(out), out, row.names = NULL)
}

# BH-FDR correction applied to the exposure row only
apply_fdr <- function(coef_table, group_vars = "sex") {
  coef_table |>
    filter(term == EXPOSURE) |>
    group_by(across(all_of(group_vars))) |>
    mutate(
      p_fdr   = p.adjust(p_value, method = "BH"),
      sig_fdr = p_fdr < 0.05
    ) |>
    ungroup()
}

# =============================================================================
# 2. Data loading and preprocessing
# =============================================================================

cat("=== Loading data ===\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE)

prepare_data <- function(data) {
  data |>
    rename(
      tyg       = TyG,
      smoking   = smoking_status,
      bmi       = BMI0,
      sbp       = SBP_mean,
      dbp       = DBP_mean,
      scr_raw   = Creatinine,
      crp       = CRP,
      hdl       = HDL.cholesterol,
      ldl       = LDL,
      chol      = Cholesterol,
      drinking  = drinking.status,
      education = Education_Level
    ) |>
    mutate(
      sex = case_when(
        sex %in% c("F", "Female", "female", "0") ~ "Female",
        sex %in% c("M", "Male",   "male",   "1") ~ "Male",
        TRUE ~ as.character(sex)
      ),
      smoking   = as.character(smoking),
      drinking  = as.character(drinking),
      education = as.character(education),
      TDI       = as.numeric(TDI),
      scr_raw   = scr_raw / 88.4,          # μmol/L → mg/dL
      egfr      = calc_egfr(scr_raw, age, sex)
    ) |>
    select(-scr_raw) |>
    select(all_of(ALL_VARS))
}

analysis_data   <- prepare_data(raw)
continuous_cols <- c("tyg", BIOMARKERS, "age", "TDI")

# Primary analysis: remove 5-SD outliers
data_main <- analysis_data |>
  mutate(across(all_of(continuous_cols), remove_outliers)) |>
  filter(complete.cases(.))

# Sensitivity analysis: retain outliers
data_sens <- analysis_data |>
  filter(complete.cases(.))

cat(sprintf("Raw data      : %d rows\n", nrow(raw)))
cat(sprintf("Primary (5-SD): %d rows\n", nrow(data_main)))
cat(sprintf("Sensitivity   : %d rows\n", nrow(data_sens)))

sex_strat      <- split(data_main, data_main[[SEX_VAR]])
sex_strat_sens <- split(data_sens, data_sens[[SEX_VAR]])

# =============================================================================
# 3. TyG–biomarker regression
# =============================================================================

fit_tyg_models <- function(sex_strats) {
  map_dfr(set_names(SEX_GROUPS), function(sg) {
    tibble(Biomarker = BIOMARKERS) |>
      mutate(mod = map(Biomarker, ~ lm(
        reformulate(c(EXPOSURE, COVARIATES), response = .x),
        data = sex_strats[[sg]], na.action = na.exclude
      )))
  }, .id = "sex")
}

get_coefs <- function(mdtb) {
  mdtb |>
    mutate(mod = map(mod, extract_coefs)) |>
    unnest(mod)
}

cat("\n=== TyG regression (primary analysis) ===\n")
mdtb          <- fit_tyg_models(sex_strat)
bmi_coefs     <- get_coefs(mdtb)
bmi_coefs_fdr <- apply_fdr(bmi_coefs)
cat("TyG effects after BH-FDR correction:\n")
print(bmi_coefs_fdr)

# Interaction tests: TyG × covariate
cat("\n=== TyG × covariate interaction tests ===\n")
INTERACTION_VARS <- c("age", "smoking", "drinking")

interaction_effects <- map_dfr(set_names(SEX_GROUPS), function(sg) {
  map_dfr(set_names(BIOMARKERS), function(bm) {
    int_terms <- paste0(EXPOSURE, ":", INTERACTION_VARS)
    mod <- lm(
      reformulate(c(EXPOSURE, COVARIATES, int_terms), response = bm),
      data = sex_strat[[sg]], na.action = na.exclude
    )
    extract_coefs(mod)
  }, .id = "Biomarker")
}, .id = "sex")

interaction_sig <- interaction_effects |>
  filter(grepl(":", term)) |>
  group_by(sex, Biomarker) |>
  mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  filter(p_fdr < 0.05)

if (nrow(interaction_sig) > 0) {
  cat("Significant interactions (FDR < 0.05):\n")
  print(interaction_sig)
} else {
  cat("No significant interactions found.\n")
}

# =============================================================================
# 4. Standardised residuals
# =============================================================================

get_residuals <- function(sex_strats, mdtb) {
  map(set_names(SEX_GROUPS), function(sg) {
    sub <- mdtb |> filter(sex == sg)
    resid_mat <- map(sub$mod, ~ scale(as.vector(resid(.x)))) |>
      set_names(sub$Biomarker) |>
      as_tibble()
    tibble(eid = sex_strats[[sg]][[ID_COL]]) |> bind_cols(resid_mat)
  })
}

cat("\n=== Computing standardised residuals ===\n")
residuals_by_sex <- get_residuals(sex_strat, mdtb)

# =============================================================================
# 5. GMM soft assignment
# =============================================================================

gmm_fit <- function(data_mat, cluster_pars) {
  mus     <- map(cluster_pars, "center")
  covmats <- map(cluster_pars, "cov")
  pdfs    <- mapply(function(mu, cov) dmvnorm(data_mat, mu, cov), mus, covmats)
  k       <- length(mus)
  ws      <- rep(1 / k, k)

  L     <- pdfs %*% diag(ws)
  probs <- L / rowSums(L)
  ll    <- sum(log(rowSums(L)))
  delta <- 1; iter <- 1

  while (delta > 1e-10 && iter <= 100) {
    ws     <- colSums(probs) / sum(probs)
    L      <- pdfs %*% diag(ws)
    probs  <- L / rowSums(L)
    ll_new <- sum(log(rowSums(L)))
    delta  <- ll_new - ll
    ll     <- ll_new
    iter   <- iter + 1
  }

  colnames(probs) <- names(mus)
  message("  GMM converged in ", iter, " iterations.")
  list(probs = as.data.frame(probs), weights = ws)
}

# =============================================================================
# 6. UMAP + Leiden clustering
# =============================================================================

run_umap_clustering <- function(resid_mat) {
  set.seed(RANDOM_SEED)

  eids <- resid_mat[[ID_COL]]
  xmat <- resid_mat[, names(resid_mat) != ID_COL]
  n    <- nrow(xmat)
  nn   <- max(10, round(10 + 15 * (log10(n) - 4)))

  message(sprintf("  n = %d, nn = %d, seed = %d", n, nn, RANDOM_SEED))

  # UMAP embedding
  umap_res <- umap(
    xmat, n_components = 2, n_neighbors = nn,
    nn_method = "annoy", n_trees = 100,
    n_sgd_threads = "auto", init = "pca",
    n_epochs = 500, binary_edge_weights = TRUE,
    dens_scale = 1, ret_extra = c("model", "fgraph"),
    verbose = FALSE
  )

  # Build nearest-neighbour graph
  g <- graph_from_adjacency_matrix(umap_res$fgraph, mode = "undirected")
  V(g)$name <- eids

  # Two-stage Leiden clustering
  set.seed(RANDOM_SEED)
  init_clus  <- cluster_leading_eigen(g)
  set.seed(RANDOM_SEED)
  final_clus <- cluster_leiden(
    g, objective_function = "modularity",
    initial_membership = init_clus$membership,
    n_iterations = 500
  )

  n_clus <- length(sizes(final_clus))
  message(sprintf("  %d clusters found; modularity = %.3f", n_clus, final_clus$quality))

  # Eigenvector-centrality-weighted Gaussian parameters per cluster
  member_list <- split(
    data.frame(
      eid     = as.numeric(names(membership(final_clus))),
      cluster = as.numeric(membership(final_clus))
    ),
    as.numeric(membership(final_clus))
  )

  cluster_pars <- map(member_list, function(cl) {
    subg   <- induced_subgraph(g, as.character(cl$eid))
    ec     <- eigen_centrality(subg, directed = FALSE)$vector
    cl_mat <- xmat[match(as.numeric(names(ec)), eids), , drop = FALSE]
    cov.wt(cl_mat, wt = as.numeric(ec))
  })

  # Add baseline concordant (BC) cluster: mean = 0, cov = I
  cluster_pars[["0"]] <- list(
    center = setNames(rep(0, ncol(xmat)), colnames(xmat)),
    cov    = diag(ncol(xmat)) |>
      `dimnames<-`(list(colnames(xmat), colnames(xmat)))
  )

  # GMM
  gmm_res <- gmm_fit(xmat, cluster_pars)
  for (i in seq_along(cluster_pars))
    cluster_pars[[i]]$weight <- gmm_res$weights[[i]]
  names(cluster_pars) <- paste0("cluster_", names(cluster_pars))

  prob_df <- data.frame(eid = eids, gmm_res$probs)
  colnames(prob_df) <- gsub("^X", "cluster_", colnames(prob_df))

  embed_df <- data.frame(eid = eids, umap_res$embedding)
  colnames(embed_df) <- c(ID_COL, "UMAP1", "UMAP2")

  list(
    probs       = prob_df,
    clusters    = cluster_pars,
    modularity  = final_clus$quality,
    embedding   = embed_df,
    random_seed = RANDOM_SEED
  )
}

cat("\n=== UMAP clustering (primary analysis, seed = 42) ===\n")
clustering_results <- map(set_names(SEX_GROUPS), function(sg) {
  message("Processing: ", sg)
  run_umap_clustering(residuals_by_sex[[sg]])
})

# =============================================================================
# 7. Cluster-level TyG effects (probability-weighted regression)
# =============================================================================

get_cluster_tyg_effects <- function(clus_res, sex_strats) {
  map_dfr(set_names(SEX_GROUPS), function(sg) {
    probs <- clus_res[[sg]]$probs |>
      pivot_longer(-eid, names_to = "cluster", values_to = "w")
    inner_join(probs, sex_strats[[sg]], by = "eid") |>
      group_by(cluster) |>
      nest() |>
      mutate(results = map(data, function(d) {
        map_dfr(set_names(BIOMARKERS), function(bm) {
          lm(reformulate(c(EXPOSURE, COVARIATES), bm),
             data = d, weights = d$w) |>
            extract_coefs()
        }, .id = "Biomarker")
      })) |>
      select(-data) |>
      unnest(results)
  }, .id = "sex")
}

cat("\n=== Cluster-level TyG effects ===\n")
cluster_tyg_effects <- get_cluster_tyg_effects(clustering_results, sex_strat)

cluster_tyg_effects_fdr <- cluster_tyg_effects |>
  filter(term == EXPOSURE) |>
  group_by(sex, cluster) |>
  mutate(p_fdr = p.adjust(p_value, method = "BH"), sig_fdr = p_fdr < 0.05) |>
  ungroup()

# =============================================================================
# 8. Sensitivity analysis (outliers retained)
# =============================================================================

cat("\n=== Sensitivity analysis (outliers retained) ===\n")

mdtb_sens          <- fit_tyg_models(sex_strat_sens)
bmi_coefs_sens     <- get_coefs(mdtb_sens)
bmi_coefs_sens_fdr <- apply_fdr(bmi_coefs_sens)

residuals_sens  <- get_residuals(sex_strat_sens, mdtb_sens)

clustering_sens <- map(set_names(SEX_GROUPS), function(sg) {
  message("Sensitivity: ", sg)
  run_umap_clustering(residuals_sens[[sg]])
})

cluster_tyg_effects_sens     <- get_cluster_tyg_effects(clustering_sens, sex_strat_sens)
cluster_tyg_effects_sens_fdr <- cluster_tyg_effects_sens |>
  filter(term == EXPOSURE) |>
  group_by(sex, cluster) |>
  mutate(p_fdr = p.adjust(p_value, method = "BH"), sig_fdr = p_fdr < 0.05) |>
  ungroup()

sensitivity_comparison <- left_join(
  bmi_coefs_fdr      |> select(sex, Biomarker, Est_main = Estimate, p_fdr_main = p_fdr),
  bmi_coefs_sens_fdr |> select(sex, Biomarker, Est_sens = Estimate, p_fdr_sens = p_fdr),
  by = c("sex", "Biomarker")
) |>
  mutate(Consistent = sign(Est_main) == sign(Est_sens))

cat("\n=== Robustness check ===\n")
cat(sprintf(
  "Direction consistency (primary vs sensitivity): %.1f%%\n",
  100 * mean(sensitivity_comparison$Consistent, na.rm = TRUE)
))

# =============================================================================
# 9. Save results
# =============================================================================

complete_results <- list(
  random_seed             = RANDOM_SEED,
  covariates_used         = COVARIATES,
  # Primary analysis
  sex_stratified          = sex_strat,
  bmi_coefs               = bmi_coefs,
  bmi_coefs_fdr           = bmi_coefs_fdr,
  interaction_effects     = interaction_effects,
  clustering_results      = clustering_results,
  cluster_tyg_effects     = cluster_tyg_effects,
  cluster_tyg_effects_fdr = cluster_tyg_effects_fdr,
  # Sensitivity analysis
  sex_stratified_sens          = sex_strat_sens,
  bmi_coefs_sens_fdr           = bmi_coefs_sens_fdr,
  clustering_results_sens      = clustering_sens,
  cluster_tyg_effects_sens_fdr = cluster_tyg_effects_sens_fdr,
  sensitivity_comparison       = sensitivity_comparison
)

save(complete_results,
     file = file.path(OUTPUT_DIR, "TyG_clustering_results.RData"))

cat("\n Done. Results saved to TyG_clustering_results.RData\n")
