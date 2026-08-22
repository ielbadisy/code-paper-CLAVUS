source("R/00_setup.R")

# ── Phase 1: Data acquisition ─────────────────────────────────────────────────
source("R/01_download_clinvar.R")
source("R/02_parse_clinvar_releases.R")
source("R/03_build_longitudinal_vus_dataset.R")
source("R/04_annotate_variants.R")

# ── Phase 2: Extended feature integration ─────────────────────────────────────
# Downloads AlphaMissense full-genome scores, gnomAD gene constraint (LOEUF/pLI),
# ClinGen gene-disease validity, and derives ACMG criterion proxy flags.
# Gracefully skips any source that is unavailable; writes warnings to logs.
source("R/03b_download_external_features.R")

# ── Phase 3: Feature engineering (picks up extended features from 03b) ────────
source("R/05_feature_engineering.R")

# ── Phase 4: Model training — expanded learner pool (13 learners + ensemble) ──
# Adds BART and KNN to the original 11-learner pool.
source("R/06_train_compare_learners.R")

# ── Phase 5: Model selection ──────────────────────────────────────────────────
source("R/07_select_best_model.R")

# ── Phase 6: Validation ───────────────────────────────────────────────────────
source("R/08_validate_best_model.R")

# ── Phase 7: Calibration, DCA, enhanced metrics ───────────────────────────────
# Reliability diagrams, ECE/MCE/ICI/HL, Decision Curve Analysis,
# DeLong pairwise tests, NRI/IDI, bootstrapped CIs, confidence-tier analysis.
source("R/08b_calibration_dca_metrics.R")

# ── Phase 8: Subgroup and equity analysis ─────────────────────────────────────
# Performance by variant type, molecular consequence, review status, gene class,
# KM-style cumulative reclassification curve by risk group.
source("R/08c_subgroup_analysis.R")

# ── Phase 9: Interpretability ─────────────────────────────────────────────────
source("R/09_interpretability.R")

# ── Phase 10: Treatment effect estimation — multiple estimators ───────────────
# G-computation (funcml), IPW with stabilized weights, AIPW (doubly-robust),
# simple one-step TMLE, sign-direction agreement check, E-values.
source("R/10_ate_estimation.R")

# ── Phase 11: Tables and Nature Medicine publication-quality figures ───────────
source("R/11_tables_figures.R")

# ── Phase 12: Render all reports ─────────────────────────────────────────────
render_qmd <- function(path) {
  if (requireNamespace("quarto", quietly = TRUE)) {
    quarto::quarto_render(path, quiet = FALSE)
  } else {
    system2("quarto", c("render", path))
  }
}

render_qmd("reports/protocol.qmd")
render_qmd("reports/analysis_report.qmd")
render_qmd("reports/manuscript_skeleton.qmd")

# Render the Nature Medicine submission manuscript (primary deliverable)
if (file.exists("reports/nature_medicine_manuscript.qmd")) {
  render_qmd("reports/nature_medicine_manuscript.qmd")
} else {
  log_message("outputs/logs/pipeline_warnings.txt",
    "nature_medicine_manuscript.qmd not found; skipping Nature Medicine PDF render.")
}

log_message("outputs/logs/pipeline_warnings.txt",
  sprintf("Full pipeline completed successfully at %s.", format(Sys.time())))
message("Pipeline complete. Key outputs in outputs/ and reports/.")
