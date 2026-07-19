# TyG–Biomarker Discordance Profiling

Unsupervised profiling of TyG–biomarker discordance in sex-stratified cohorts using UMAP, Leiden clustering, and Gaussian Mixture Models (GMM). Downstream analyses cover cluster-specific TyG effects, prevalent disease associations, medication patterns, and longitudinal MACE/T2D risk.

## Pipeline Overview

| Script | Description | Key Output |
|--------|-------------|------------|
| `01_clustering.R` | Data prep · UMAP · Leiden · GMM soft assignment | `TyG_clustering_results.RData` |
| `02_biomarker_profiles.R` | Cluster weight bar charts + residual scatter plots | `TyG_cluster_overview.png/pdf` |
| `03_overall_tyg_effects.R` | Forest plots of TyG → biomarker effects (overall) | `TyG_overall_effects.png/pdf` |
| `04_cluster_tyg_effects.R` | Cluster-specific weighted regression + Excel table | `TyG_cluster_effects.png/pdf`, `Supplementary_Table_Cluster_TyG_Effects.xlsx` |
| `05_prevalent_disease.R` | ALR-logistic regression for prevalent diseases | `TyG_prevalent_disease_forest.png/pdf` |
| `06_medication.R` | ALR-logistic regression for medication use | `TyG_medication_forest.png/pdf` |
| `07_survival.R` | ALR-Cox models for MACE and incident T2D | `TyG_survival_forest.png/pdf` |

## Requirements

```r
install.packages(c(
  "dplyr", "tibble", "tidyr", "purrr", "readr",
  "uwot", "igraph", "mvtnorm",
  "ggplot2", "cowplot", "patchwork", "ggh4x", "scales",
  "openxlsx", "survival"
))
```

## Data

Place `imp_data1.csv` in the working directory before running any script. The file must contain the columns listed in `01_clustering.R` (`ALL_VARS`, disease indicators, medication flags, and survival date fields).

## Usage

Run scripts in order; each script loads the `.RData` produced by `01_clustering.R`.

```r
source("01_clustering.R")   # ~10–20 min depending on sample size
source("02_biomarker_profiles.R")
source("03_overall_tyg_effects.R")
source("04_cluster_tyg_effects.R")
source("05_prevalent_disease.R")
source("06_medication.R")
source("07_survival.R")
```

All outputs are written to the working directory (set `OUTPUT_DIR` in each script to redirect).

## Cluster Labels

| Label | Phenotype |
|-------|-----------|
| BC | Baseline Concordant (reference) |
| DHT | Discordant Hypertensive |
| DRI | Discordant Renal Insufficiency |
| DHL | Discordant High HDL-Lipid |
| DOB | Discordant Obesity |
| DPL | Discordant Pro-Atherogenic Lipid |
| DIS | Discordant Inflammatory |
