# CLAVUS -- Step 2: Train, select, validate, interpret, and estimate adjusted effects
#
# Trains nine supervised learners (funcml) under time-aware rolling
# cross-validation, selects the gradient boosting model by temporal-test
# AUPRC, evaluates it under temporal and gene-held-out validation, computes
# calibration/decision-curve/subgroup diagnostics, runs interpretability
# analyses (VIP, SHAP, PDP/ICE/ALE), and estimates adjusted treatment effects
# (IPW/AIPW/TMLE) for four clinically interpretable exposures.
#
# Run from the project root, after 01_build_dataset.R (or using the curated
# dataset already provided in data/curated/):
#   Rscript 02_run_ml_pipeline.R

source("R/00_setup.R")
source("R/06_train_compare_learners.R")
source("R/07_select_best_model.R")
source("R/08_validate_best_model.R")
source("R/08b_calibration_dca_metrics.R")
source("R/08c_subgroup_analysis.R")
source("R/09_interpretability.R")
source("R/10_ate_estimation.R")
