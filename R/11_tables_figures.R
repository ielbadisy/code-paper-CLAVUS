## 11_tables_figures.R
## Nature Medicine publication-quality tables and figures
## Replaces the previous version; keeps all table logic, adds polished figures.

source("R/00_setup.R")
load_or_stop("ggplot2"); library(ggplot2)
load_or_stop("patchwork"); library(patchwork)

## ── Colorblind-safe Okabe-Ito palette ────────────────────────────────────────
OI <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#000000")
names(OI) <- c("orange", "sky", "green", "yellow",
               "blue", "vermillion", "pink", "black")

## ── Helper: save figure as PNG + PDF ─────────────────────────────────────────
save_fig <- function(p, stem, width_px = 1200, height_px = 900) {
  png_path <- file.path("outputs/figures", paste0(stem, ".png"))
  pdf_path <- file.path("outputs/figures", paste0(stem, ".pdf"))
  png(png_path, width = width_px, height = height_px, res = 120)
  print(p)
  dev.off()
  pdf(pdf_path, width = width_px / 120, height = height_px / 120)
  print(p)
  dev.off()
  invisible(NULL)
}

## ── Nature Medicine theme ─────────────────────────────────────────────────────
nm_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text        = ggplot2::element_text(size = 10, color = "black"),
      axis.title       = ggplot2::element_text(size = 12),
      legend.text      = ggplot2::element_text(size = 10),
      legend.title     = ggplot2::element_text(size = 10, face = "bold"),
      plot.title       = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "#444444"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = 10, face = "bold"),
      panel.border     = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5)
    )
}

## ═══════════════════════════════════════════════════════════════════════════════
## SECTION 1 — TABLES (preserve original logic verbatim)
## ═══════════════════════════════════════════════════════════════════════════════

dat        <- readRDS("data/curated/vus_resolution_dataset.rds")
perf       <- read.csv("outputs/tables/model_performance.csv")
validation <- read.csv("outputs/tables/best_model_validation.csv")
vi         <- read.csv("outputs/tables/variable_importance_best_model.csv")
ate        <- read.csv("outputs/tables/ate_results.csv")
pred_all   <- read.csv("outputs/predictions/temporal_test_predictions_all_learners.csv")
validity   <- read.csv("outputs/tables/dataset_validity_checks.csv")
best_name  <- readLines("outputs/models/best_model_name.txt")

## Table 1: cohort characteristics
table1_vars <- c("time_to_reclassification_days", "baseline_number_submitters",
                 "baseline_number_conditions", "days_since_first_seen",
                 "cadd_score", "revel_score", "gnomad_af")
table1_vars <- table1_vars[vapply(dat[, table1_vars, drop = FALSE],
                                  function(x) !all(is.na(x)), logical(1))]
table1 <- aggregate(
  stats::as.formula(paste("cbind(", paste(table1_vars, collapse = ","), ") ~ final_class_group")),
  dat, function(x) round(mean(x, na.rm = TRUE), 3)
)
overall_values <- as.list(vapply(dat[, names(table1)[-1]],
                                  function(x) round(mean(x, na.rm = TRUE), 3), numeric(1)))
overall <- data.frame(c(list(final_class_group = "Overall"), overall_values), check.names = FALSE)
names(overall) <- names(table1)
utils::write.csv(rbind(overall, table1),
                 "outputs/tables/table1_cohort_characteristics.csv", row.names = FALSE)

## Table 2: missingness summary
missingness <- data.frame(
  variable  = names(dat),
  n_missing = vapply(dat, function(x) sum(is.na(x)), integer(1)),
  pct_missing = round(vapply(dat, function(x) mean(is.na(x)), numeric(1)) * 100, 2)
)
utils::write.csv(missingness, "outputs/tables/table2_missingness_summary.csv", row.names = FALSE)

## Tables 3–8
utils::write.csv(perf,       "outputs/tables/table3_performance_comparison.csv",     row.names = FALSE)
utils::write.csv(validation, "outputs/tables/table4_best_model_validation.csv",      row.names = FALSE)
utils::write.csv(head(vi, 30), "outputs/tables/table5_top_predictors.csv",           row.names = FALSE)
utils::write.csv(ate,        "outputs/tables/table6_ate_adjusted_effects.csv",       row.names = FALSE)
utils::write.csv(validity,   "outputs/tables/table7_dataset_validity_checks.csv",    row.names = FALSE)
if (file.exists("outputs/tables/sensitivity_subgroup_validation.csv")) {
  sensitivity <- read.csv("outputs/tables/sensitivity_subgroup_validation.csv")
  utils::write.csv(sensitivity, "outputs/tables/table8_sensitivity_subgroup_validation.csv",
                   row.names = FALSE)
}

## Leakage audit
recipe  <- readRDS("data/curated/feature_recipe.rds")
forbidden <- c("final_class", "final_class_group", "final_date", "final_review_status",
               "final_release_month", "label", "time_to_reclassification_days",
               "best_model_score", "risk_group")
leakage_audit <- data.frame(
  check  = c("no_outcome_or_final_fields_in_predictors",
             "external_predictors_have_real_coverage_or_missingness_indicator",
             "temporal_test_after_train",
             "same_variant_not_in_train_and_test"),
  passed = c(
    length(intersect(recipe$predictor_names, forbidden)) == 0,
    {
      ext <- intersect(recipe$predictor_names,
                       c("cadd_score", "revel_score", "gnomad_af",
                         "alphamissense_score", "spliceai_score"))
      all(vapply(ext, function(nm) sum(!is.na(dat[[nm]])) >= 50, logical(1))) &&
        all(c("is_missing_gnomad_af", "is_missing_cadd", "is_missing_revel",
               "is_missing_alphamissense", "is_missing_spliceai") %in% recipe$predictor_names)
    },
    all(validity$passed[validity$check == "temporal_split_order_valid"]),
    all(validity$passed[validity$check == "train_test_variant_disjoint"])
  )
)
utils::write.csv(leakage_audit, "outputs/tables/leakage_audit.csv", row.names = FALSE)

## Calibration diagnostics
calibration_summary <- perf[, c("learner", "Brier_score", "calibration_intercept",
                                 "calibration_slope", "ECE")]
calibration_summary$selected_best_model <- calibration_summary$learner == best_name
utils::write.csv(calibration_summary, "outputs/tables/calibration_diagnostics.csv",
                 row.names = FALSE)

best_cal <- calibration_summary[calibration_summary$selected_best_model, ]
validity_status <- if (all(validity$passed)) {
  "All dataset validity checks passed before modeling."
} else {
  "At least one dataset validity check failed."
}
slope_note <- if (best_cal$calibration_slope < 0.8) {
  "Calibration slope is below 1, suggesting over-dispersed risk scores on the temporal test split."
} else if (best_cal$calibration_slope > 1.25) {
  "Calibration slope is above 1, suggesting under-dispersed risk scores on the temporal test split."
} else {
  "Calibration slope is close to 1, suggesting acceptable spread of predicted probabilities on the temporal test split."
}
ece_note <- if (best_cal$ECE <= 0.05) {
  "ECE is low; binned observed risks are close to predicted probabilities."
} else if (best_cal$ECE <= 0.10) {
  "ECE is moderate; calibration is usable for prioritization but should still be inspected visually."
} else {
  "ECE is high; the score should be interpreted cautiously and recalibrated before operational use."
}
ate_primary <- ate[ate$exposure == "High computational deleteriousness evidence" &
                     ate$method == "funcml_g_computation", ]
ate_gene    <- ate[ate$exposure == "Gene pathogenic-enrichment above median at t0" &
                     ate$method == "funcml_g_computation", ]
ate_note <- if (nrow(ate_primary) > 0 && !is.na(ate_primary$estimate)) {
  paste("The high-computational-evidence adjusted association is",
        round(ate_primary$estimate, 4),
        "on the risk-difference scale using real external annotation coverage where available.")
} else if (nrow(ate_gene) > 0 && !is.na(ate_gene$estimate)) {
  paste("The gene-enrichment adjusted association is",
        round(ate_gene$estimate, 4),
        "on the risk-difference scale and is based on real ClinVar-derived variables.")
} else {
  "The adjusted association analysis was infeasible on available real data."
}
interpretation <- data.frame(
  topic = c("dataset_validity", "best_model", "calibration_slope",
            "calibration_error", "adjusted_effects", "interpretability_scope"),
  interpretation = c(
    validity_status,
    paste("Selected best model:", best_name,
          "based on temporal-test AUPRC with Brier/calibration tie-breakers."),
    slope_note, ece_note, ate_note,
    "Variable-importance and diagnostic plots explain prediction behavior, not causal pathogenic mechanisms or definitive clinical truth."
  )
)
utils::write.csv(interpretation, "outputs/tables/result_interpretation.csv", row.names = FALSE)

## Tuning plots (funcml)
tuning_path <- "outputs/models/funcml_tuning_results.rds"
if (file.exists(tuning_path)) {
  tuning <- readRDS(tuning_path)
  for (nm in names(tuning)) {
    png(file.path("outputs/figures", paste0("funcml_tuning_", nm, ".png")),
        width = 1000, height = 700)
    print(plot(tuning[[nm]]))
    dev.off()
  }
}

writeLines(capture.output(sessionInfo()), "outputs/logs/session_info.txt")

## ═══════════════════════════════════════════════════════════════════════════════
## SECTION 2 — NATURE MEDICINE FIGURES
## ═══════════════════════════════════════════════════════════════════════════════

cat("\n[11_tables_figures] Generating Nature Medicine figures ...\n")

## ── FIGURE 1: Study design + Cohort flow ─────────────────────────────────────

## Panel A: schematic timeline using annotate / geom_segment
boxes <- data.frame(
  x    = c(0.05, 0.40, 0.75),
  y    = c(0.50, 0.50, 0.50),
  xmin = c(0.00, 0.35, 0.70),
  xmax = c(0.20, 0.55, 0.90),
  ymin = c(0.35, 0.35, 0.35),
  ymax = c(0.65, 0.65, 0.65),
  label = c("2020-01\nBaseline VUS", "2021-2026\nFollow-up releases",
            "VUS reclassification\noutcome")
)

pA <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = boxes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#D0E4F7", color = OI["blue"], linewidth = 0.8
  ) +
  ggplot2::geom_text(
    data = boxes,
    ggplot2::aes(x = x, y = y, label = label),
    size = 3.2, lineheight = 1.2
  ) +
  ggplot2::annotate(
    "segment",
    x = 0.20, xend = 0.35, y = 0.50, yend = 0.50,
    arrow = grid::arrow(length = grid::unit(0.025, "npc"), type = "closed"),
    color = OI["blue"], linewidth = 0.9
  ) +
  ggplot2::annotate(
    "segment",
    x = 0.55, xend = 0.70, y = 0.50, yend = 0.50,
    arrow = grid::arrow(length = grid::unit(0.025, "npc"), type = "closed"),
    color = OI["blue"], linewidth = 0.9
  ) +
  ggplot2::annotate(
    "text", x = 0.275, y = 0.60,
    label = "Predictor\nextraction", size = 2.8, color = "#333333"
  ) +
  ggplot2::annotate(
    "text", x = 0.625, y = 0.60,
    label = "Outcome\nascertainment", size = 2.8, color = "#333333"
  ) +
  ggplot2::annotate(
    "text", x = 0.45, y = 0.15,
    label = "Endpoint: ClinVar reclassification to P/LP or B/LB",
    size = 3.0, color = "#555555", fontface = "italic"
  ) +
  ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
  ggplot2::labs(title = "Study timeline") +
  nm_theme() +
  ggplot2::theme(
    axis.text  = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.line  = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

## Panel B: cohort flow bar chart with real numbers
n_baseline  <- nrow(dat)
n_resolved  <- sum(!is.na(dat$final_class_group))
n_analytic  <- sum(!is.na(dat$label))
## all variants in follow-up (we have outcome for all, so n_followup = n_baseline)
n_followup  <- n_baseline

flow_df <- data.frame(
  step = factor(
    c("Baseline VUS (2020-01)",
      "Variants in follow-up releases",
      "Resolved to P/LP or B/LB",
      "Final analytic cohort"),
    levels = c("Final analytic cohort",
               "Resolved to P/LP or B/LB",
               "Variants in follow-up releases",
               "Baseline VUS (2020-01)")
  ),
  n = c(n_baseline, n_followup, n_resolved, n_analytic)
)

pB <- ggplot2::ggplot(flow_df, ggplot2::aes(x = n, y = step)) +
  ggplot2::geom_col(fill = OI["blue"], width = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = formatC(n, format = "d", big.mark = ",")),
    hjust = -0.08, size = 3.5
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) formatC(x, format = "d", big.mark = ","),
    expand = ggplot2::expansion(mult = c(0, 0.18))
  ) +
  ggplot2::labs(x = "Number of variants", y = NULL, title = "Cohort flow") +
  nm_theme()

fig1 <- pA + pB +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"))
  )

save_fig(fig1, "figure1_study_design_cohort", width_px = 1800, height_px = 750)
cat("[fig1] Study design + cohort flow saved.\n")

## ── FIGURE 2: Learner comparison ─────────────────────────────────────────────

perf_ord <- perf[order(perf$AUPRC, decreasing = FALSE), ]
perf_ord$learner_f <- factor(perf_ord$learner, levels = perf_ord$learner)
perf_ord$is_best   <- perf_ord$learner == best_name
fill_vals <- c("TRUE" = OI["vermillion"], "FALSE" = "#AABBCC")

pA2 <- ggplot2::ggplot(perf_ord, ggplot2::aes(x = AUPRC, y = learner_f, fill = is_best)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(x = "AUPRC", y = NULL) +
  nm_theme()

pB2 <- ggplot2::ggplot(perf_ord, ggplot2::aes(x = AUROC, y = learner_f, fill = is_best)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(x = "AUROC", y = NULL) +
  nm_theme() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

pC2 <- ggplot2::ggplot(perf_ord, ggplot2::aes(x = Brier_score, y = learner_f, fill = is_best)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
  ggplot2::scale_x_reverse(limits = c(0.30, 0), breaks = seq(0, 0.30, 0.05)) +
  ggplot2::labs(x = "Brier score (lower = better)", y = NULL,
                caption = paste0("Best model: ", best_name, " (highlighted)")) +
  nm_theme() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

fig2 <- pA2 + pB2 + pC2 +
  patchwork::plot_annotation(
    tag_levels = "a",
    title      = "Learner comparison on temporal test set",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.tag   = ggplot2::element_text(size = 14, face = "bold")
    )
  )

save_fig(fig2, "figure2_learner_comparison")
cat("[fig2] Learner comparison saved.\n")

## ── FIGURE 3: ROC and PR curves ──────────────────────────────────────────────

## Identify score columns (exclude label/id/raw columns)
score_cols <- setdiff(
  names(pred_all)[!grepl("_raw$", names(pred_all))],
  c("variant_id", "label", "final_class_group")
)

compute_curves <- function(score, y, learner) {
  thresholds <- sort(unique(c(seq(0, 1, length.out = 201), score)))
  do.call(rbind, lapply(thresholds, function(th) {
    pr <- as.integer(score >= th)
    tp <- sum(pr == 1 & y == 1)
    fp <- sum(pr == 1 & y == 0)
    fn <- sum(pr == 0 & y == 1)
    tn <- sum(pr == 0 & y == 0)
    data.frame(
      learner   = learner,
      threshold = th,
      tpr       = tp / max(tp + fn, 1),
      fpr       = fp / max(fp + tn, 1),
      precision = tp / max(tp + fp, 1),
      recall    = tp / max(tp + fn, 1)
    )
  }))
}

curve_list <- lapply(score_cols, function(nm) {
  compute_curves(pred_all[[nm]], pred_all$label, nm)
})
curve_points <- do.call(rbind, curve_list)

## AUROC per learner for legend annotation
auroc_vec <- vapply(score_cols, function(nm) {
  perf_row <- perf[perf$learner == nm, ]
  if (nrow(perf_row) == 1) perf_row$AUROC else NA_real_
}, numeric(1))

legend_labels <- setNames(
  paste0(score_cols, sprintf("  (AUROC=%.3f)", auroc_vec)),
  score_cols
)

## Color: best model in vermillion, others in shades of blue/grey
n_learners    <- length(score_cols)
other_colors  <- grDevices::colorRampPalette(c("#A0C4D8", "#4A7DA8"))(n_learners - 1)
all_colors    <- setNames(
  c(other_colors, OI["vermillion"]),
  c(setdiff(score_cols, best_name), best_name)
)
linewidths    <- setNames(
  c(rep(0.6, n_learners - 1), 1.4),
  c(setdiff(score_cols, best_name), best_name)
)

pA3 <- ggplot2::ggplot(curve_points,
    ggplot2::aes(x = fpr, y = tpr, color = learner,
                 linewidth = learner, group = learner)) +
  ggplot2::geom_line() +
  ggplot2::geom_abline(slope = 1, intercept = 0,
                       linetype = "dashed", color = "gray60", linewidth = 0.4) +
  ggplot2::scale_color_manual(values = all_colors, labels = legend_labels,
                              name = NULL) +
  ggplot2::scale_linewidth_manual(values = linewidths, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1), name = "False positive rate") +
  ggplot2::scale_y_continuous(limits = c(0, 1), name = "True positive rate") +
  ggplot2::labs(title = "ROC curves") +
  nm_theme() +
  ggplot2::theme(legend.position = "bottom",
                 legend.text = ggplot2::element_text(size = 8)) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 2))

pB3 <- ggplot2::ggplot(curve_points,
    ggplot2::aes(x = recall, y = precision, color = learner,
                 linewidth = learner, group = learner)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = all_colors, guide = "none") +
  ggplot2::scale_linewidth_manual(values = linewidths, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1), name = "Recall") +
  ggplot2::scale_y_continuous(limits = c(0, 1), name = "Precision") +
  ggplot2::labs(title = "Precision-recall curves") +
  nm_theme()

fig3 <- pA3 + pB3 +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"))
  )

save_fig(fig3, "figure3_roc_pr_curves", width_px = 1800, height_px = 850)
cat("[fig3] ROC/PR curves saved.\n")

## ── FIGURE 4: Calibration of best model ──────────────────────────────────────

best_pred_path <- "outputs/predictions/best_model_temporal_test_predictions.csv"

if (file.exists(best_pred_path)) {

  best_pred <- read.csv(best_pred_path)
  score_col <- if ("best_model_score" %in% names(best_pred)) "best_model_score" else
    setdiff(names(best_pred), c("variant_id", "label", "final_class_group", "risk_group"))[1]
  score  <- best_pred[[score_col]]
  y      <- best_pred$label

  ## 10-bin reliability diagram
  n_bins <- 10
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bin_id <- findInterval(score, breaks, rightmost.closed = TRUE)
  bin_id <- pmin(bin_id, n_bins)

  cal_df <- do.call(rbind, lapply(seq_len(n_bins), function(b) {
    idx <- bin_id == b
    if (sum(idx) == 0) return(NULL)
    data.frame(
      bin      = b,
      mid      = (breaks[b] + breaks[b + 1]) / 2,
      mean_pred = mean(score[idx]),
      obs_freq  = mean(y[idx]),
      n        = sum(idx)
    )
  }))

  ## ECE
  ece_val <- if ("ECE" %in% names(best_cal) && !is.na(best_cal$ECE)) {
    round(best_cal$ECE, 4)
  } else {
    round(sum(cal_df$n * abs(cal_df$mean_pred - cal_df$obs_freq)) / sum(cal_df$n), 4)
  }

  pA4 <- ggplot2::ggplot(cal_df, ggplot2::aes(x = mean_pred, y = obs_freq)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "gray50", linewidth = 0.7) +
    ggplot2::geom_smooth(
      data = data.frame(score = score, y = y),
      ggplot2::aes(x = score, y = y),
      inherit.aes = FALSE,
      method = "loess", formula = y ~ x, se = FALSE,
      color = OI["green"], linewidth = 0.9, linetype = "solid"
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n), color = OI["blue"], alpha = 0.85) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "n in bin") +
    ggplot2::annotate(
      "text", x = 0.78, y = 0.10,
      label = paste0("ECE = ", ece_val),
      size = 3.8, fontface = "italic", color = "#333333"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), name = "Mean predicted probability") +
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "Observed event frequency") +
    ggplot2::labs(title = "Reliability diagram",
                  subtitle = paste0("Best model: ", best_name, "  |  10 bins")) +
    nm_theme()

  ## Panel B: score density by class
  dens_df <- data.frame(
    score = score,
    class = ifelse(y == 1, "P/LP", "B/LB")
  )
  pB4 <- ggplot2::ggplot(dens_df, ggplot2::aes(x = score, fill = class, color = class)) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.8) +
    ggplot2::scale_fill_manual(
      values = c("P/LP" = OI["vermillion"], "B/LB" = OI["sky"]),
      name   = "Final class"
    ) +
    ggplot2::scale_color_manual(
      values = c("P/LP" = OI["vermillion"], "B/LB" = OI["sky"]),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), name = "Predicted probability (P/LP)") +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::labs(
      title    = "Score distribution by outcome class",
      subtitle = sprintf("N = %d  (P/LP = %d, B/LB = %d)",
                         nrow(best_pred), sum(y == 1), sum(y == 0))
    ) +
    nm_theme() +
    ggplot2::theme(legend.position = c(0.82, 0.82))

  fig4 <- pA4 + pB4 +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"))
    )

  save_fig(fig4, "figure4_calibration")
  cat("[fig4] Calibration saved.\n")

} else {
  message("[fig4] SKIP: best_model_temporal_test_predictions.csv not found.")
}

## ── FIGURE 5: Variable importance ────────────────────────────────────────────

vi_path <- "outputs/tables/variable_importance_best_model.csv"

if (file.exists(vi_path)) {

  vi_df <- read.csv(vi_path)

  ## Tidy variable names for display
  tidy_varname <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("is ", "", x)
    x <- gsub("baseline ", "baseline\n", x)
    x <- gsub("gene ", "gene\n", x)
    tools::toTitleCase(trimws(x))
  }

  top20 <- head(vi_df[order(vi_df$importance, decreasing = TRUE), ], 20)
  top20$var_label <- factor(tidy_varname(top20$variable),
                            levels = rev(tidy_varname(top20$variable)))

  pA5 <- ggplot2::ggplot(top20, ggplot2::aes(x = importance, y = var_label)) +
    ggplot2::geom_col(fill = OI["blue"], width = 0.65) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)),
                                 name = "Permutation importance") +
    ggplot2::labs(y = NULL, title = "Top 20 predictors",
                  subtitle = paste0("Model: ", best_name)) +
    nm_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

  shap_path <- "outputs/figures/funcml_interpret_shap_best_model.png"

  if (file.exists(shap_path)) {
    ## Embed SHAP image in a ggplot panel using annotation_raster or grid
    shap_img <- tryCatch(
      png::readPNG(shap_path),
      error = function(e) NULL
    )
    if (!is.null(shap_img)) {
      pB5 <- ggplot2::ggplot() +
        ggplot2::annotation_raster(shap_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "SHAP values (best model)") +
        nm_theme() +
        ggplot2::theme(
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.line  = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank()
        )
    } else {
      shap_img <- NULL
    }
  } else {
    shap_img <- NULL
  }

  if (is.null(shap_img)) {
    ## Simplified direction-of-effect panel: top 10 features, signed by correlation with score
    if (file.exists(best_pred_path)) {
      bp    <- read.csv(best_pred_path)
      dat2  <- readRDS("data/curated/vus_resolution_dataset.rds")
      merged_vi <- merge(bp[, c("variant_id", "best_model_score")],
                         dat2, by = "variant_id", all.x = TRUE)
      top10_vars <- as.character(head(vi_df[order(vi_df$importance, decreasing = TRUE), "variable"], 10))
      top10_vars <- intersect(top10_vars, names(merged_vi))
      dirs <- vapply(top10_vars, function(v) {
        x  <- as.numeric(merged_vi[[v]])
        ok <- !is.na(x) & !is.na(merged_vi$best_model_score)
        if (sum(ok) < 10) return(0)
        cor(x[ok], merged_vi$best_model_score[ok], method = "spearman")
      }, numeric(1))
      dir_df <- data.frame(
        variable  = factor(tidy_varname(top10_vars), levels = rev(tidy_varname(top10_vars))),
        direction = dirs
      )
      pB5 <- ggplot2::ggplot(dir_df, ggplot2::aes(x = direction, y = variable,
                                                    fill = direction > 0)) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::scale_fill_manual(
          values = c("TRUE" = OI["vermillion"], "FALSE" = OI["sky"]),
          labels = c("Higher = lower P/LP risk", "Higher = higher P/LP risk"),
          name   = "Effect direction"
        ) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
        ggplot2::labs(x = "Spearman correlation with score", y = NULL,
                      title = "Feature direction of effect (top 10)") +
        nm_theme() +
        ggplot2::theme(legend.position = "bottom",
                       axis.text.y = ggplot2::element_text(size = 9))
    } else {
      pB5 <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "SHAP figure not available", size = 5) +
        ggplot2::theme_void()
    }
  }

  fig5 <- pA5 + pB5 +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"))
    )

  save_fig(fig5, "figure5_variable_importance")
  cat("[fig5] Variable importance saved.\n")

} else {
  message("[fig5] SKIP: variable_importance_best_model.csv not found.")
}

## ── FIGURE 6: Decision curve analysis ────────────────────────────────────────

dca_path <- "outputs/tables/dca_results.csv"

if (file.exists(dca_path)) {

  dca_wide <- read.csv(dca_path)
  # Reshape from wide (model / treat_all / treat_none) to long format
  dca_long_from_wide <- function(d) {
    cols <- c("model", "treat_all", "treat_none")
    labels <- c("Model", "Treat all", "Treat none")
    present <- cols[cols %in% names(d)]
    do.call(rbind, lapply(seq_along(present), function(i) {
      data.frame(threshold = d$threshold, net_benefit = d[[present[i]]],
                 strategy = labels[i], stringsAsFactors = FALSE)
    }))
  }

  make_dca_panel <- function(dca_sub, outcome_label) {
    long <- if ("strategy" %in% names(dca_sub)) dca_sub else dca_long_from_wide(dca_sub)
    strategy_colors <- c("Model" = OI["blue"], "Treat all" = OI["vermillion"], "Treat none" = "gray60")
    strategy_lty    <- c("Model" = "solid", "Treat all" = "dashed", "Treat none" = "dotted")
    present_strats  <- unique(long$strategy)
    ggplot2::ggplot(long,
      ggplot2::aes(x = threshold, y = net_benefit, color = strategy, linetype = strategy)) +
      ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", linewidth = 0.4) +
      ggplot2::scale_x_continuous(name = "Probability threshold",
                                   labels = function(x) paste0(round(x * 100), "%")) +
      ggplot2::scale_y_continuous(name = "Net benefit") +
      ggplot2::scale_color_manual(values = strategy_colors[present_strats], name = "Strategy") +
      ggplot2::scale_linetype_manual(values = strategy_lty[present_strats], name = "Strategy") +
      ggplot2::labs(title = outcome_label) +
      nm_theme() +
      ggplot2::theme(legend.position = "bottom")
  }

  has_plp <- "outcome" %in% names(dca_wide) && any(grepl("P_LP|P/LP|plp", dca_wide$outcome, ignore.case = TRUE))
  has_blb <- "outcome" %in% names(dca_wide) && any(grepl("B_LB|B/LB|blb", dca_wide$outcome, ignore.case = TRUE))

  if (has_plp && has_blb) {
    dca_plp <- dca_wide[grepl("P_LP|P/LP|plp", dca_wide$outcome, ignore.case = TRUE), ]
    dca_blb <- dca_wide[grepl("B_LB|B/LB|blb", dca_wide$outcome, ignore.case = TRUE), ]
    pA6 <- make_dca_panel(dca_plp, "DCA: P/LP prediction")
    pB6 <- make_dca_panel(dca_blb, "DCA: B/LB prediction")
  } else {
    pA6 <- make_dca_panel(dca_wide, "Decision curve analysis")
    pB6 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "B/LB DCA panel not available", size = 4.5, color = "gray40") +
      ggplot2::theme_void()
  }

  fig6 <- pA6 + pB6 +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = "Decision curve analysis",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.tag   = ggplot2::element_text(size = 14, face = "bold")
      )
    )

  save_fig(fig6, "figure6_decision_curve")
  cat("[fig6] Decision curve analysis saved.\n")

} else {
  message("[fig6] SKIP: dca_results.csv not found — creating placeholder.")
  fig6 <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5,
                      label = "Figure 6: Decision curve analysis\n(run DCA module first to generate dca_results.csv)",
                      size = 5, color = "gray40", hjust = 0.5) +
    ggplot2::theme_void()
  save_fig(fig6, "figure6_decision_curve")
}

## ── FIGURE 7: ATE forest plot ─────────────────────────────────────────────────

ate_path <- "outputs/tables/ate_results.csv"

if (file.exists(ate_path)) {

  ate_full <- read.csv(ate_path)

  ## Determine sign agreement across methods per exposure
  ate_agree <- tapply(sign(ate_full$estimate), ate_full$exposure, function(s) {
    length(unique(s[!is.na(s)])) <= 1
  })

  ate_full$sign_agree <- ate_agree[ate_full$exposure]
  ate_full$exposure_short <- substring(ate_full$exposure, 1, 40)

  ## Sort exposures by mean estimate
  exp_order <- names(sort(tapply(ate_full$estimate, ate_full$exposure, mean, na.rm = TRUE)))
  ate_full$exposure_f <- factor(ate_full$exposure_short,
                                 levels = substring(exp_order, 1, 40))

  method_colors <- setNames(
    OI[seq_len(length(unique(ate_full$method)))],
    unique(ate_full$method)
  )

  fig7 <- ggplot2::ggplot(ate_full,
    ggplot2::aes(x = estimate, y = exposure_f,
                 xmin = conf_low, xmax = conf_high,
                 color = method, shape = method)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray50", linewidth = 0.7) +
    ggplot2::geom_errorbarh(height = 0.3, linewidth = 0.8,
                             position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_point(size = 3,
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_color_manual(values = method_colors, name = "Method") +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 8)[seq_len(length(unique(ate_full$method)))],
                                 name = "Method") +
    ggplot2::scale_x_continuous(name = "Adjusted risk difference (95% CI)") +
    ggplot2::labs(
      y        = NULL,
      title    = "Average treatment effect estimates",
      subtitle = "Risk-difference scale; G-computation adjustment"
    ) +
    ggplot2::facet_grid(sign_agree ~ ., labeller = ggplot2::as_labeller(
      c("TRUE" = "Sign agreement across methods",
        "FALSE" = "Inconsistent sign across methods")
    ), scales = "free_y", space = "free_y") +
    nm_theme() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(size = 9))

  save_fig(fig7, "figure7_ate_forest_plot", width_px = 1300, height_px = 900)
  cat("[fig7] ATE forest plot saved.\n")

} else {
  message("[fig7] SKIP: ate_results.csv not found.")
}

## ── FIGURE 8: Subgroup forest plot + KM curves ───────────────────────────────

sg_path <- "outputs/tables/subgroup_performance.csv"
km_path <- "outputs/figures/km_by_risk_group.png"

if (file.exists(sg_path)) {

  sg_df <- read.csv(sg_path)
  sg_df <- sg_df[!is.na(sg_df$AUROC), ]
  sg_df$label_n <- paste0(sg_df$subgroup, "  (n=", sg_df$n, ")")

  ## Sort: Overall at top, rest by AUROC ascending
  is_overall  <- sg_df$subgroup == "Overall"
  sg_other    <- sg_df[!is_overall, ]
  sg_over     <- sg_df[is_overall, ]
  sg_other    <- sg_other[order(sg_other$AUROC, decreasing = FALSE), ]
  sg_sorted   <- rbind(sg_other, sg_over)
  sg_sorted$label_n <- factor(sg_sorted$label_n, levels = sg_sorted$label_n)

  overall_auroc_sg <- sg_sorted$AUROC[sg_sorted$subgroup == "Overall"]

  pA8 <- ggplot2::ggplot(
    sg_sorted,
    ggplot2::aes(x = AUROC, y = label_n,
                 xmin = AUROC_lower, xmax = AUROC_upper,
                 color = subgroup == "Overall")
  ) +
    ggplot2::geom_vline(xintercept = overall_auroc_sg,
                        linetype = "dashed", color = "#555555", linewidth = 0.7) +
    ggplot2::geom_errorbarh(height = 0.3, linewidth = 0.75, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = OI["blue"], "TRUE" = OI["vermillion"]),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(0.45, 1.0),
                                 breaks = seq(0.5, 1.0, by = 0.1)) +
    ggplot2::labs(
      x        = "AUROC (95% bootstrap CI)",
      y        = NULL,
      title    = "Subgroup AUROC",
      subtitle = paste0("Overall = ", round(overall_auroc_sg, 3),
                        " (dashed line)")
    ) +
    nm_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  ## Panel B: load KM figure or recreate if exists
  if (file.exists(km_path) && requireNamespace("png", quietly = TRUE)) {
    km_img <- tryCatch(png::readPNG(km_path), error = function(e) NULL)
    if (!is.null(km_img)) {
      pB8 <- ggplot2::ggplot() +
        ggplot2::annotation_raster(km_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "Cumulative reclassification by risk group") +
        nm_theme() +
        ggplot2::theme(
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.line  = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank()
        )
    } else {
      km_img <- NULL
    }
  } else {
    km_img <- NULL
  }

  if (is.null(km_img)) {
    pB8 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "KM figure: run 08c_subgroup_analysis.R first",
                        size = 4, color = "gray40") +
      ggplot2::theme_void()
  }

  fig8 <- pA8 + pB8 +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = "Subgroup performance and clinical actionability",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.tag   = ggplot2::element_text(size = 14, face = "bold")
      )
    )

  save_fig(fig8, "figure8_subgroup_km", width_px = 1800, height_px = 900)
  cat("[fig8] Subgroup forest + KM saved.\n")

} else {
  message("[fig8] SKIP: subgroup_performance.csv not found — run 08c_subgroup_analysis.R first.")
}

## ── Done ─────────────────────────────────────────────────────────────────────
cat("\n[11_tables_figures] All tables and figures complete.\n")
