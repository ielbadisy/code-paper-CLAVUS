## 08c_subgroup_analysis.R
## Subgroup and equity analysis on the temporal test set
## Outputs: subgroup_performance.csv, subgroup_auroc_forest_plot.png, km_by_risk_group.png

source("R/00_setup.R")
load_or_stop("ggplot2"); library(ggplot2)

## ── 0. Data loading ─────────────────────────────────────────────────────────

dat  <- readRDS("data/curated/vus_resolution_dataset.rds")
pred <- read.csv("outputs/predictions/best_model_temporal_test_predictions.csv",
                 stringsAsFactors = FALSE)

## Merge predictions with dataset to obtain subgroup variables
keep_cols <- c("variant_id", "molecular_consequence", "variant_type",
               "baseline_review_status", "gene_symbol",
               "gene_pathogenic_fraction_at_t0", "time_to_reclassification_days",
               "is_missense", "is_splice_region", "is_nonsense", "is_frameshift")
keep_cols <- intersect(keep_cols, names(dat))

df <- merge(pred, dat[, keep_cols, drop = FALSE], by = "variant_id", all.x = TRUE)
df$score  <- df$best_model_score
df$y      <- df$label  # 1 = P/LP, 0 = B/LB

cat(sprintf("[subgroup] Merged dataset: %d rows, %d columns\n", nrow(df), ncol(df)))

## ── 1. Metric helper functions ───────────────────────────────────────────────

compute_auroc <- function(score, y) {
  if (length(unique(y)) < 2) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r   <- rank(score)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

compute_auprc <- function(score, y) {
  if (length(unique(y)) < 2) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y_s <- y[ord]
  tp  <- cumsum(y_s)
  fp  <- cumsum(1 - y_s)
  prec <- tp / (tp + fp)
  rec  <- tp / sum(y)
  ## trapezoidal AUC under PR curve
  delta_rec <- c(rec[1], diff(rec))
  sum(prec * delta_rec)
}

compute_brier <- function(score, y) mean((score - y)^2)

compute_sens_spec <- function(score, y, threshold = 0.5) {
  pred_class <- as.integer(score >= threshold)
  tp <- sum(pred_class == 1 & y == 1)
  tn <- sum(pred_class == 0 & y == 0)
  fp <- sum(pred_class == 1 & y == 0)
  fn <- sum(pred_class == 0 & y == 1)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  c(sensitivity = sens, specificity = spec)
}

## Bootstrap AUROC with 95% CI
bootstrap_auroc <- function(score, y, B = 200, seed = 42) {
  set.seed(seed)
  n    <- length(y)
  boot <- vapply(seq_len(B), function(i) {
    idx <- sample(n, n, replace = TRUE)
    compute_auroc(score[idx], y[idx])
  }, numeric(1))
  boot <- boot[!is.na(boot)]
  if (length(boot) < 10) return(c(auroc = compute_auroc(score, y), ci_low = NA, ci_high = NA))
  c(auroc   = compute_auroc(score, y),
    ci_low  = quantile(boot, 0.025),
    ci_high = quantile(boot, 0.975))
}

## Wrapper: compute all metrics for a subset
subgroup_metrics <- function(score, y, label, n_total, B = 200) {
  if (n_total < 50 || length(unique(y)) < 2) return(NULL)
  boot  <- bootstrap_auroc(score, y, B = B)
  ss    <- compute_sens_spec(score, y)
  data.frame(
    subgroup    = label,
    n           = n_total,
    n_P_LP      = sum(y == 1),
    n_B_LB      = sum(y == 0),
    AUROC       = round(boot["auroc"],   4),
    AUROC_lower = round(boot["ci_low"],  4),
    AUROC_upper = round(boot["ci_high"], 4),
    AUPRC       = round(compute_auprc(score, y), 4),
    Brier       = round(compute_brier(score, y), 4),
    sensitivity = round(ss["sensitivity"], 4),
    specificity = round(ss["specificity"], 4),
    row.names   = NULL
  )
}

## ── 2. Overall ───────────────────────────────────────────────────────────────

results <- list()

results[["Overall"]] <- subgroup_metrics(
  df$score, df$y,
  label   = "Overall",
  n_total = nrow(df),
  B       = 200
)

## ── 3. Subgroup: Molecular consequence ──────────────────────────────────────
## Derive a 5-level variable from the boolean flags and molecular_consequence

derive_mol_cons <- function(d) {
  out <- rep("other", nrow(d))
  if ("is_missense"     %in% names(d)) out[d$is_missense     == 1] <- "missense"
  if ("is_frameshift"   %in% names(d)) out[d$is_frameshift   == 1] <- "frameshift"
  if ("is_nonsense"     %in% names(d)) out[d$is_nonsense     == 1] <- "nonsense"
  if ("is_splice_region" %in% names(d)) out[d$is_splice_region == 1] <- "splice_region"
  out
}

df$mol_cons_group <- derive_mol_cons(df)

mol_levels <- c("missense", "splice_region", "frameshift", "nonsense", "other")
for (lv in mol_levels) {
  idx <- df$mol_cons_group == lv
  sub <- df[idx, ]
  key <- paste0("mol_cons_", lv)
  results[[key]] <- subgroup_metrics(
    sub$score, sub$y,
    label   = paste0("Mol. consequence: ", lv),
    n_total = nrow(sub),
    B       = 200
  )
}

## ── 4. Subgroup: Variant type ────────────────────────────────────────────────

collapse_vtype <- function(vt) {
  snv_pat <- "single nucleotide"
  del_pat <- "^deletion$"
  ins_pat <- "^insertion$"
  dup_pat <- "^duplication$"
  out <- rep("other", length(vt))
  out[grepl(snv_pat, vt, ignore.case = TRUE)] <- "single nucleotide variant"
  out[grepl(del_pat, vt, ignore.case = TRUE, perl = TRUE)] <- "deletion"
  out[grepl(ins_pat, vt, ignore.case = TRUE, perl = TRUE)] <- "insertion"
  out[grepl(dup_pat, vt, ignore.case = TRUE, perl = TRUE)] <- "duplication"
  out
}

df$vtype_group <- collapse_vtype(df$variant_type)

vtype_levels <- c("single nucleotide variant", "deletion", "insertion", "duplication", "other")
for (lv in vtype_levels) {
  idx <- df$vtype_group == lv
  sub <- df[idx, ]
  key <- paste0("vtype_", gsub(" ", "_", lv))
  results[[key]] <- subgroup_metrics(
    sub$score, sub$y,
    label   = paste0("Variant type: ", lv),
    n_total = nrow(sub),
    B       = 200
  )
}

## ── 5. Subgroup: Baseline review status ─────────────────────────────────────

collapse_review <- function(rs) {
  out <- rep("other", length(rs))
  out[grepl("no assertion", rs, ignore.case = TRUE)]           <- "no assertion criteria"
  out[grepl("single submitter", rs, ignore.case = TRUE)]       <- "criteria: single submitter"
  out[grepl("multiple submitters", rs, ignore.case = TRUE)]    <- "criteria: multiple submitters"
  out[grepl("expert panel", rs, ignore.case = TRUE)]           <- "expert panel"
  out
}

df$review_group <- collapse_review(df$baseline_review_status)

review_levels <- c("no assertion criteria", "criteria: single submitter",
                   "criteria: multiple submitters", "expert panel")
for (lv in review_levels) {
  idx <- df$review_group == lv
  sub <- df[idx, ]
  key <- paste0("review_", gsub("[: /]", "_", lv))
  results[[key]] <- subgroup_metrics(
    sub$score, sub$y,
    label   = paste0("Review status: ", lv),
    n_total = nrow(sub),
    B       = 200
  )
}

## ── 6. Subgroup: Gene class ──────────────────────────────────────────────────

hereditary_cancer_panel <- c(
  "BRCA1", "BRCA2", "TP53", "MLH1", "MSH2", "MSH6", "PMS2",
  "APC", "MUTYH", "PTEN", "CDH1", "STK11", "PALB2", "ATM",
  "CHEK2", "BRIP1", "RAD51C", "RAD51D"
)

classify_gene <- function(gs, gpf) {
  out <- rep("other genes", length(gs))
  out[gpf >= 0.3 & !gs %in% hereditary_cancer_panel] <- "high pathogenic fraction (>=0.3)"
  out[gs %in% hereditary_cancer_panel] <- "well-known disease genes"
  out
}

df$gene_class <- classify_gene(
  df$gene_symbol,
  ifelse(is.na(df$gene_pathogenic_fraction_at_t0), 0, df$gene_pathogenic_fraction_at_t0)
)

gene_levels <- c("well-known disease genes", "high pathogenic fraction (>=0.3)", "other genes")
for (lv in gene_levels) {
  idx <- df$gene_class == lv
  sub <- df[idx, ]
  key <- paste0("gene_class_", gsub("[^a-zA-Z0-9]", "_", lv))
  results[[key]] <- subgroup_metrics(
    sub$score, sub$y,
    label   = paste0("Gene class: ", lv),
    n_total = nrow(sub),
    B       = 200
  )
}

## ── 7. Compile subgroup table ────────────────────────────────────────────────

subgroup_perf <- do.call(rbind, Filter(Negate(is.null), results))
rownames(subgroup_perf) <- NULL

utils::write.csv(subgroup_perf, "outputs/tables/subgroup_performance.csv", row.names = FALSE)
cat(sprintf("[subgroup] Subgroup performance table saved (%d rows).\n", nrow(subgroup_perf)))

## ── 8. KM-style cumulative reclassification by risk group ───────────────────

## Merge time_to_reclassification_days from dataset
if ("time_to_reclassification_days" %in% names(df) &&
    !all(is.na(df$time_to_reclassification_days))) {

  km_df <- df[!is.na(df$time_to_reclassification_days) &
                !is.na(df$risk_group), ]

  ## Colorblind-safe Okabe-Ito colors
  okabe_ito <- c(
    "likely B/LB resolution"     = "#56B4E9",
    "remain uncertain / abstain" = "#E69F00",
    "likely P/LP resolution"     = "#D55E00"
  )

  ## Empirical CDF per risk group
  compute_ecdf_df <- function(times, group_name) {
    t_sorted <- sort(times)
    n        <- length(t_sorted)
    cum_frac <- seq_len(n) / n
    data.frame(
      time     = t_sorted,
      cum_frac = cum_frac,
      risk_group = group_name
    )
  }

  groups     <- unique(km_df$risk_group)
  ecdf_parts <- lapply(groups, function(g) {
    sub_times <- km_df$time_to_reclassification_days[km_df$risk_group == g]
    if (length(sub_times) == 0) return(NULL)
    ## Prepend origin point
    origin <- data.frame(time = 0, cum_frac = 0, risk_group = g)
    rbind(origin, compute_ecdf_df(sub_times, g))
  })
  ecdf_df <- do.call(rbind, Filter(Negate(is.null), ecdf_parts))

  ## Group sample sizes for legend
  n_per_group <- table(km_df$risk_group)
  ecdf_df$risk_group_label <- paste0(
    ecdf_df$risk_group, " (n=",
    n_per_group[ecdf_df$risk_group], ")"
  )

  p_km <- ggplot2::ggplot(ecdf_df,
    ggplot2::aes(x = time, y = cum_frac,
                 color = risk_group, group = risk_group)) +
    ggplot2::geom_step(linewidth = 1.1) +
    ggplot2::scale_color_manual(
      values = okabe_ito,
      labels = paste0(names(n_per_group), " (n=", as.integer(n_per_group), ")"),
      name   = "Predicted risk group"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(round(x / 365.25, 1), " yr"),
      name   = "Days to reclassification"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      name   = "Cumulative fraction reclassified"
    ) +
    ggplot2::labs(
      title    = "Cumulative reclassification by predicted risk group",
      subtitle = sprintf("Temporal test set  (N = %d)", nrow(km_df))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 13, face = "bold")
    )

  ## Use scales if available, otherwise suppress y-axis label formatting
  if (!requireNamespace("scales", quietly = TRUE)) {
    p_km <- p_km +
      ggplot2::scale_y_continuous(limits = c(0, 1), name = "Cumulative fraction reclassified")
  }

  png("outputs/figures/km_by_risk_group.png", width = 1000, height = 750, res = 120)
  print(p_km)
  dev.off()
  cat("[subgroup] KM cumulative reclassification figure saved.\n")

} else {
  cat("[subgroup] WARNING: time_to_reclassification_days not available in merged data; skipping KM figure.\n")
}

## ── 9. Forest plot of subgroup AUROC ────────────────────────────────────────

forest_df <- subgroup_perf[!is.na(subgroup_perf$AUROC), ]
forest_df$label_n <- paste0(forest_df$subgroup, "  (n=", forest_df$n, ")")

## Sort by AUROC descending, Overall at top
overall_idx  <- forest_df$subgroup == "Overall"
forest_other <- forest_df[!overall_idx, ]
forest_over  <- forest_df[overall_idx, ]
forest_other <- forest_other[order(forest_other$AUROC, decreasing = FALSE), ]
forest_df    <- rbind(forest_other, forest_over)
forest_df$label_n <- factor(forest_df$label_n, levels = forest_df$label_n)

overall_auroc <- forest_df$AUROC[forest_df$subgroup == "Overall"]

p_forest <- ggplot2::ggplot(
  forest_df,
  ggplot2::aes(x = AUROC, y = label_n,
               xmin = AUROC_lower, xmax = AUROC_upper,
               color = subgroup == "Overall")
) +
  ggplot2::geom_vline(xintercept = overall_auroc, linetype = "dashed",
                      color = "#555555", linewidth = 0.7) +
  ggplot2::geom_errorbarh(height = 0.25, linewidth = 0.8, na.rm = TRUE) +
  ggplot2::geom_point(size = 3, na.rm = TRUE) +
  ggplot2::scale_color_manual(
    values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
    guide  = "none"
  ) +
  ggplot2::scale_x_continuous(limits = c(0.45, 1.0),
                               breaks = seq(0.5, 1.0, by = 0.1)) +
  ggplot2::labs(
    x        = "AUROC (95% bootstrap CI, B=200)",
    y        = NULL,
    title    = "Subgroup AUROC — temporal test set",
    subtitle = paste0("Vertical dashed line = overall AUROC (",
                      round(overall_auroc, 3), ")")
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_text(size = 9),
    axis.text.x  = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(size = 12),
    plot.title   = ggplot2::element_text(size = 13, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 10, color = "#444444")
  )

png("outputs/figures/subgroup_auroc_forest_plot.png",
    width = 1100, height = 900, res = 120)
print(p_forest)
dev.off()
cat("[subgroup] Subgroup AUROC forest plot saved.\n")

## ── Done ─────────────────────────────────────────────────────────────────────
cat("\n[08c_subgroup_analysis] All subgroup analyses complete.\n")
cat("  - outputs/tables/subgroup_performance.csv\n")
cat("  - outputs/figures/km_by_risk_group.png\n")
cat("  - outputs/figures/subgroup_auroc_forest_plot.png\n")
