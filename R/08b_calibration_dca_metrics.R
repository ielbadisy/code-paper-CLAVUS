## =============================================================================
## 08b_calibration_dca_metrics.R
## Advanced calibration, decision curve analysis, and enhanced statistical
## metrics. Runs AFTER 08_validate_best_model.R.
## =============================================================================

source("R/00_setup.R")

## ── packages ─────────────────────────────────────────────────────────────────
load_or_stop("ggplot2");  library(ggplot2)
load_or_stop("patchwork"); library(patchwork)

## ── colour palette (Okabe-Ito, colorblind-safe) ──────────────────────────────
OI_PAL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#000000")

## ── learner names ─────────────────────────────────────────────────────────────
LEARNERS <- c("xgboost", "bart", "ranger", "gbm", "earth",
              "glmnet", "svm", "nnet", "knn")

LEARNER_LABELS <- c(
  xgboost = "XGBoost",
  bart    = "BART",
  ranger  = "Random Forest",
  gbm     = "GBM",
  earth   = "MARS",
  glmnet  = "Elastic Net",
  svm     = "Radial SVM",
  nnet    = "Neural Network",
  knn     = "k-NN"
)

## ── load predictions ──────────────────────────────────────────────────────────
preds_all <- utils::read.csv(
  "outputs/predictions/temporal_test_predictions_all_learners.csv",
  stringsAsFactors = FALSE
)

best_name <- tryCatch(
  trimws(readLines("outputs/models/best_model_name.txt")[1]),
  error = function(e) "gradient_boosting"
)

## outcome vector (binary 0/1)
y_all <- preds_all$label

## ── helpers ──────────────────────────────────────────────────────────────────

## guard: return NA if predictions are degenerate
is_degenerate <- function(p) {
  length(unique(round(p, 6))) < 2 || all(is.na(p))
}

## AUROC (Wilcoxon-Mann-Whitney)
auroc <- function(y, p) {
  if (is_degenerate(p) || length(unique(y)) < 2) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## AUPRC
auprc <- function(y, p) {
  if (is_degenerate(p) || length(unique(y)) < 2) return(NA_real_)
  ord <- order(p, decreasing = TRUE)
  yy  <- y[ord]
  tp  <- cumsum(yy == 1)
  fp  <- cumsum(yy == 0)
  prec <- tp / pmax(tp + fp, 1)
  rec  <- tp / max(sum(yy == 1), 1)
  sum(diff(c(0, rec)) * prec)
}

## Brier score
brier <- function(y, p) mean((y - p)^2, na.rm = TRUE)

## logit
logit <- function(p) log(p / (1 - p))

## safe logit (clip first)
safe_logit <- function(p) logit(pmin(pmax(p, 1e-9), 1 - 1e-9))

## ECE (equal-width bins)
ece <- function(y, p, n_bins = 10) {
  if (is_degenerate(p)) return(NA_real_)
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bins   <- cut(p, breaks = breaks, include.lowest = TRUE)
  g <- split(seq_along(y), bins)
  n <- length(y)
  sum(vapply(g, function(idx) {
    if (length(idx) == 0) return(0)
    abs(mean(p[idx]) - mean(y[idx])) * length(idx) / n
  }, numeric(1)), na.rm = TRUE)
}

## MCE (maximum calibration error across bins)
mce <- function(y, p, n_bins = 10) {
  if (is_degenerate(p)) return(NA_real_)
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bins   <- cut(p, breaks = breaks, include.lowest = TRUE)
  g <- split(seq_along(y), bins)
  max(vapply(g, function(idx) {
    if (length(idx) == 0) return(0)
    abs(mean(p[idx]) - mean(y[idx]))
  }, numeric(1)), na.rm = TRUE)
}

## Calibration-in-the-large
citl <- function(y, p) mean(p, na.rm = TRUE) - mean(y, na.rm = TRUE)

## Calibration slope
calib_slope <- function(y, p) {
  if (is_degenerate(p) || length(unique(y)) < 2) return(NA_real_)
  lp <- safe_logit(p)
  tryCatch(coef(glm(y ~ lp, family = binomial))[["lp"]], error = function(e) NA_real_)
}

## Hosmer-Lemeshow test (10 groups, decile-based)
hosmer_lemeshow <- function(y, p, g = 10) {
  if (is_degenerate(p) || length(unique(y)) < 2)
    return(list(statistic = NA_real_, p_value = NA_real_))
  ord    <- order(p)
  y_ord  <- y[ord]
  p_ord  <- p[ord]
  groups <- cut(seq_along(p_ord), breaks = g, labels = FALSE)
  hl_stat <- sum(vapply(seq_len(g), function(k) {
    idx <- which(groups == k)
    if (length(idx) == 0) return(0)
    obs1 <- sum(y_ord[idx] == 1)
    obs0 <- length(idx) - obs1
    exp1 <- sum(p_ord[idx])
    exp0 <- length(idx) - exp1
    v1 <- if (exp1 > 0) (obs1 - exp1)^2 / exp1 else 0
    v0 <- if (exp0 > 0) (obs0 - exp0)^2 / exp0 else 0
    v1 + v0
  }, numeric(1)))
  df <- g - 2
  list(statistic = hl_stat, p_value = pchisq(hl_stat, df = df, lower.tail = FALSE))
}

## ICI (Integrated Calibration Index via Loess)
ici_metric <- function(y, p, ...) {
  if (is_degenerate(p) || length(unique(y)) < 2) return(NA_real_)
  tryCatch({
    fit  <- loess(y ~ p, span = 0.75, family = "symmetric", ...)
    pred <- predict(fit, newdata = data.frame(p = p))
    mean(abs(pred - p), na.rm = TRUE)
  }, error = function(e) NA_real_)
}

## ── loess calibration curve data ─────────────────────────────────────────────
loess_calib_df <- function(y, p, n_grid = 200) {
  if (is_degenerate(p) || length(unique(y)) < 2) return(NULL)
  tryCatch({
    fit   <- loess(y ~ p, span = 0.75, family = "symmetric")
    pgrid <- seq(min(p), max(p), length.out = n_grid)
    pred  <- predict(fit, newdata = data.frame(p = pgrid), se = TRUE)
    data.frame(
      x    = pgrid,
      y    = pmax(0, pmin(1, pred$fit)),
      ylo  = pmax(0, pmin(1, pred$fit - 1.96 * pred$se.fit)),
      yhi  = pmax(0, pmin(1, pred$fit + 1.96 * pred$se.fit))
    )
  }, error = function(e) NULL)
}

## =============================================================================
## 1. ENHANCED CALIBRATION METRICS — ALL LEARNERS
## =============================================================================

calib_rows <- lapply(LEARNERS, function(ln) {
  if (!ln %in% names(preds_all)) return(NULL)
  p <- preds_all[[ln]]
  y <- y_all
  ok <- !is.na(p) & !is.na(y)
  p  <- p[ok]; y <- y[ok]

  if (is_degenerate(p) || length(unique(y)) < 2) {
    message("  [SKIP calib] learner ", ln, ": degenerate predictions.")
    return(data.frame(
      learner      = ln,
      label        = LEARNER_LABELS[ln],
      CITL         = NA_real_, calib_slope  = NA_real_,
      ECE          = NA_real_, MCE          = NA_real_,
      HL_stat      = NA_real_, HL_pvalue    = NA_real_,
      ICI          = NA_real_
    ))
  }

  hl <- hosmer_lemeshow(y, p)
  data.frame(
    learner      = ln,
    label        = LEARNER_LABELS[ln],
    CITL         = citl(y, p),
    calib_slope  = calib_slope(y, p),
    ECE          = ece(y, p),
    MCE          = mce(y, p),
    HL_stat      = hl$statistic,
    HL_pvalue    = hl$p_value,
    ICI          = ici_metric(y, p)
  )
})

calib_df <- do.call(rbind, calib_rows[!vapply(calib_rows, is.null, logical(1))])
utils::write.csv(calib_df, "outputs/tables/calibration_metrics_all_learners.csv",
                 row.names = FALSE)
message("[08b] Calibration metrics saved.")

## =============================================================================
## 2. RELIABILITY DIAGRAMS
## =============================================================================

## ── helper: single reliability diagram gg object ─────────────────────────────
reliability_gg <- function(y, p, title = "", ece_val = NULL, ici_val = NULL,
                            base_size = 13) {
  curve_df <- loess_calib_df(y, p)
  rug_df   <- data.frame(x = p, cls = factor(y, levels = c(1, 0),
                                              labels = c("P/LP", "B/LB")))

  ann <- ""
  if (!is.null(ece_val) && !is.na(ece_val))
    ann <- paste0(ann, sprintf("ECE = %.3f", ece_val))
  if (!is.null(ici_val) && !is.na(ici_val))
    ann <- paste0(ann, if (nchar(ann) > 0) "\n" else "",
                  sprintf("ICI = %.3f", ici_val))

  g <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey50", linewidth = 0.6) +
    theme_classic(base_size = base_size) +
    scale_colour_manual(values = OI_PAL[c(1, 2)], name = "Class") +
    scale_fill_manual(values = OI_PAL[c(1, 2)], name = "Class") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Predicted probability", y = "Observed fraction",
         title = title) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = base_size, face = "bold"))

  if (!is.null(curve_df)) {
    g <- g +
      geom_ribbon(data = curve_df,
                  aes(x = x, ymin = ylo, ymax = yhi),
                  fill = OI_PAL[5], alpha = 0.25, inherit.aes = FALSE) +
      geom_line(data = curve_df,
                aes(x = x, y = y),
                colour = OI_PAL[5], linewidth = 1, inherit.aes = FALSE)
  }

  g <- g +
    geom_rug(data = rug_df, aes(x = x, colour = cls),
             sides = "b", alpha = 0.35, length = unit(0.03, "npc"),
             inherit.aes = FALSE)

  if (nchar(ann) > 0) {
    g <- g + annotate("text", x = 0.97, y = 0.05, label = ann,
                      hjust = 1, vjust = 0, size = base_size * 0.3,
                      fontface = "italic")
  }
  g
}

## ── best model: before vs after calibration ───────────────────────────────────
p_cal <- preds_all[[best_name]]
p_raw_col <- paste0(best_name, "_raw")
p_raw <- if (p_raw_col %in% names(preds_all)) preds_all[[p_raw_col]] else p_cal

ok_idx <- !is.na(p_cal) & !is.na(y_all)
y_bm   <- y_all[ok_idx]
p_bm   <- p_cal[ok_idx]
p_bm_r <- p_raw[ok_idx]

ece_after  <- ece(y_bm, p_bm)
ici_after  <- ici_metric(y_bm, p_bm)
ece_before <- ece(y_bm, p_bm_r)
ici_before <- ici_metric(y_bm, p_bm_r)

g_before <- reliability_gg(y_bm, p_bm_r,
                            title = paste0(LEARNER_LABELS[best_name],
                                           " – Before Calibration"),
                            ece_val = ece_before, ici_val = ici_before)
g_after  <- reliability_gg(y_bm, p_bm,
                            title = paste0(LEARNER_LABELS[best_name],
                                           " – After Calibration"),
                            ece_val = ece_after, ici_val = ici_after)

combined_best <- g_before + g_after + patchwork::plot_layout(ncol = 2)

png("outputs/figures/reliability_diagram_best_model.png",
    width = 1200, height = 600, res = 150)
print(combined_best)
dev.off()
message("[08b] reliability_diagram_best_model.png saved.")

## ── all-learners small multiples ─────────────────────────────────────────────
n_learners   <- length(LEARNERS)
n_cols       <- 4L
n_rows_grid  <- ceiling(n_learners / n_cols)

plot_list <- lapply(LEARNERS, function(ln) {
  if (!ln %in% names(preds_all)) return(ggplot() + theme_void())
  p  <- preds_all[[ln]]
  y  <- y_all
  ok <- !is.na(p) & !is.na(y)
  p  <- p[ok]; y <- y[ok]
  ece_v <- ece(y, p)
  ici_v <- ici_metric(y, p)
  reliability_gg(y, p, title = LEARNER_LABELS[ln],
                 ece_val = ece_v, ici_val = ici_v, base_size = 11)
})

all_grid <- patchwork::wrap_plots(plot_list, ncol = n_cols) +
  patchwork::plot_annotation(
    title   = "Reliability Diagrams – All Learners (Temporal Test Set)",
    caption = "Dashed line: perfect calibration | Shaded band: 95% Loess CI",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

png("outputs/figures/reliability_diagrams_all_learners.png",
    width = 1400, height = 900, res = 150)
print(all_grid)
dev.off()
message("[08b] reliability_diagrams_all_learners.png saved.")

## =============================================================================
## 3. DECISION CURVE ANALYSIS (from scratch)
## =============================================================================

## net benefit at threshold t for a given score vector
compute_nb <- function(y, p, thresholds) {
  n <- length(y)
  prev <- mean(y)
  vapply(thresholds, function(t) {
    pred_pos <- p >= t
    tp <- sum(pred_pos & y == 1)
    fp <- sum(pred_pos & y == 0)
    (tp / n) - (fp / n) * (t / (1 - t))
  }, numeric(1))
}

compute_nb_all <- function(y, thresholds) {
  n    <- length(y)
  prev <- mean(y)
  vapply(thresholds, function(t) {
    prev - (1 - prev) * (t / (1 - t))
  }, numeric(1))
}

thresholds <- seq(0.05, 0.95, by = 0.01)
ok_dca     <- !is.na(p_bm) & !is.na(y_bm)
y_dca      <- y_bm[ok_dca]
p_dca      <- p_bm[ok_dca]

## -- Panel A: P/LP prediction (score > t → predict P/LP) ---------------------
nb_model_plp   <- compute_nb(y_dca, p_dca, thresholds)
nb_all_plp     <- compute_nb_all(y_dca, thresholds)
nb_none_plp    <- rep(0, length(thresholds))

dca_plp <- data.frame(
  threshold   = thresholds,
  model       = nb_model_plp,
  treat_all   = nb_all_plp,
  treat_none  = nb_none_plp
)

## -- Panel B: B/LB prediction (score < (1-t) → predict B/LB) ----------------
## Flip labels so that B/LB = 1
y_blb          <- 1L - y_dca
p_blb          <- 1 - p_dca

nb_model_blb   <- compute_nb(y_blb, p_blb, thresholds)
nb_all_blb     <- compute_nb_all(y_blb, thresholds)
nb_none_blb    <- rep(0, length(thresholds))

dca_blb <- data.frame(
  threshold   = thresholds,
  model       = nb_model_blb,
  treat_all   = nb_all_blb,
  treat_none  = nb_none_blb
)

## save DCA results
dca_combined <- rbind(
  cbind(outcome = "P_LP", dca_plp),
  cbind(outcome = "B_LB", dca_blb)
)
utils::write.csv(dca_combined, "outputs/tables/dca_results.csv", row.names = FALSE)

## -- DCA plotting helper -------------------------------------------------------
dca_gg <- function(dca_df, title = "", base_size = 13) {
  long <- data.frame(
    threshold = rep(dca_df$threshold, 3),
    NB        = c(dca_df$model, dca_df$treat_all, dca_df$treat_none),
    Strategy  = rep(c("Best model", "Treat all", "Treat none"),
                    each = nrow(dca_df))
  )
  long$NB[long$NB < -0.05] <- NA   # clip extreme negatives for visual clarity
  long$Strategy <- factor(long$Strategy,
                          levels = c("Best model", "Treat all", "Treat none"))

  ## region where model > both baselines
  benefit_region <- dca_df[
    !is.na(dca_df$model) &
      dca_df$model > pmax(dca_df$treat_all, 0, na.rm = TRUE), ]

  g <- ggplot(long, aes(x = threshold, y = NB, colour = Strategy,
                         linetype = Strategy)) +
    theme_classic(base_size = base_size) +
    scale_colour_manual(
      values = c("Best model" = OI_PAL[5],
                 "Treat all"  = OI_PAL[6],
                 "Treat none" = OI_PAL[7]),
      name = "Strategy"
    ) +
    scale_linetype_manual(
      values = c("Best model" = "solid",
                 "Treat all"  = "dashed",
                 "Treat none" = "dotted"),
      name = "Strategy"
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    coord_cartesian(ylim = c(-0.05, NA)) +
    labs(x = "Threshold probability",
         y = "Net benefit",
         title = title) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = base_size, face = "bold"))

  if (nrow(benefit_region) > 0) {
    g <- g +
      annotate("rect",
               xmin = min(benefit_region$threshold),
               xmax = max(benefit_region$threshold),
               ymin = -Inf, ymax = Inf,
               alpha = 0.07, fill = OI_PAL[3]) +
      annotate("text",
               x = mean(range(benefit_region$threshold)),
               y = max(benefit_region$model, na.rm = TRUE) * 0.6,
               label = "Model advantage",
               colour = OI_PAL[3], size = 3.2, fontface = "italic")
  }
  g
}

g_dca_plp <- dca_gg(dca_plp, title = "DCA – P/LP Resolution Prediction")
g_dca_blb <- dca_gg(dca_blb, title = "DCA – B/LB Resolution Prediction")

dca_combined_plot <- g_dca_plp + g_dca_blb +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title    = paste0("Decision Curve Analysis: ", LEARNER_LABELS[best_name]),
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

png("outputs/figures/decision_curve_analysis.png",
    width = 1200, height = 600, res = 150)
print(dca_combined_plot)
dev.off()
message("[08b] decision_curve_analysis.png saved.")

## =============================================================================
## 4. DELONG PAIRWISE AUROC COMPARISON
## =============================================================================
## Implements DeLong et al. (1988) structural components approach.

## Structural components for DeLong's variance estimator
delong_components <- function(y, p) {
  y1 <- p[y == 1]; y0 <- p[y == 0]
  n1 <- length(y1); n0 <- length(y0)
  V10 <- vapply(y1, function(xi) mean((xi > y0) + 0.5 * (xi == y0)), numeric(1))
  V01 <- vapply(y0, function(xj) mean((y1 > xj) + 0.5 * (y1 == xj)), numeric(1))
  list(V10 = V10, V01 = V01, n1 = n1, n0 = n0,
       theta = mean(V10))   # = AUROC estimate
}

## DeLong test between two correlated AUROCs (same subjects)
delong_test <- function(y, p1, p2) {
  if (is_degenerate(p1) || is_degenerate(p2) || length(unique(y)) < 2)
    return(list(stat = NA_real_, pval = NA_real_))

  c1 <- delong_components(y, p1)
  c2 <- delong_components(y, p2)

  n1 <- c1$n1; n0 <- c1$n0

  ## 2x2 covariance matrix of (theta1, theta2)
  s10_11 <- var(c1$V10)
  s01_11 <- var(c1$V01)
  s10_22 <- var(c2$V10)
  s01_22 <- var(c2$V01)
  s10_12 <- cov(c1$V10, c2$V10)
  s01_12 <- cov(c1$V01, c2$V01)

  var1  <- s10_11 / n1 + s01_11 / n0
  var2  <- s10_22 / n1 + s01_22 / n0
  cov12 <- s10_12 / n1 + s01_12 / n0

  se_diff <- sqrt(var1 + var2 - 2 * cov12)
  if (se_diff <= 0) return(list(stat = NA_real_, pval = NA_real_))

  z    <- (c1$theta - c2$theta) / se_diff
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  list(stat = z, pval = pval)
}

## build pairwise matrix
valid_learners <- LEARNERS[LEARNERS %in% names(preds_all)]
n_pair        <- length(valid_learners)
n_pairs_total <- choose(n_pair, 2)

if (n_pair < 2) {
  message("[08b] Fewer than 2 valid learners for DeLong test — skipping.")
} else {

pval_mat  <- matrix(NA_real_, n_pair, n_pair,
                    dimnames = list(valid_learners, valid_learners))
stat_mat  <- pval_mat

for (i in seq_len(n_pair - 1)) {
  for (j in (i + 1):n_pair) {
    ln_i <- valid_learners[i]; ln_j <- valid_learners[j]
    p_i  <- preds_all[[ln_i]]; p_j  <- preds_all[[ln_j]]
    ok   <- !is.na(p_i) & !is.na(p_j) & !is.na(y_all)
    res  <- delong_test(y_all[ok], p_i[ok], p_j[ok])
    pval_mat[i, j] <- pval_mat[j, i] <- res$pval
    stat_mat[i, j] <- stat_mat[j, i] <- res$stat
  }
}
diag(pval_mat) <- 1

## Bonferroni correction
pval_bonf <- pmin(pval_mat * n_pairs_total, 1)

## format: p-value with asterisk if significant after correction
fmt_cell <- function(p_raw, p_b) {
  if (is.na(p_raw)) return("—")
  sig <- if (!is.na(p_b) && p_b < 0.05) "*" else ""
  paste0(sprintf("%.4f", p_raw), sig)
}

cell_mat <- matrix("—", n_pair, n_pair,
                   dimnames = list(valid_learners, valid_learners))
for (i in seq_len(n_pair)) {
  for (j in seq_len(n_pair)) {
    if (i == j) { cell_mat[i, j] <- "—"; next }
    cell_mat[i, j] <- fmt_cell(pval_mat[i, j], pval_bonf[i, j])
  }
}
diag(cell_mat) <- "—"

## reshape to long for CSV (easier to read than a wide matrix)
delong_long <- do.call(rbind, lapply(seq_len(n_pair), function(i) {
  lapply(seq_len(n_pair), function(j) {
    if (i >= j) return(NULL)
    data.frame(
      learner_A       = valid_learners[i],
      learner_B       = valid_learners[j],
      AUROC_A         = auroc(y_all[!is.na(preds_all[[valid_learners[i]]])],
                              preds_all[[valid_learners[i]]][!is.na(preds_all[[valid_learners[i]]])]),
      AUROC_B         = auroc(y_all[!is.na(preds_all[[valid_learners[j]]])],
                              preds_all[[valid_learners[j]]][!is.na(preds_all[[valid_learners[j]]])]),
      z_stat          = stat_mat[i, j],
      pvalue_raw      = pval_mat[i, j],
      pvalue_bonferroni = pval_bonf[i, j],
      significant_bonf  = !is.na(pval_bonf[i, j]) && pval_bonf[i, j] < 0.05
    )
  })
}))
delong_long <- do.call(rbind, unlist(delong_long, recursive = FALSE))

utils::write.csv(delong_long, "outputs/tables/delong_pairwise_tests.csv",
                 row.names = FALSE)
message("[08b] DeLong pairwise tests saved.")
} # end n_pair >= 2 guard

## =============================================================================
## 5. BOOTSTRAPPED CONFIDENCE INTERVALS — BEST MODEL
## =============================================================================

set.seed(20260427)
B <- 1000L

sens_at_half <- function(y, p) {
  pred <- as.integer(p >= 0.5)
  if (sum(y == 1) == 0) return(NA_real_)
  sum(pred == 1 & y == 1) / sum(y == 1)
}
spec_at_half <- function(y, p) {
  pred <- as.integer(p >= 0.5)
  if (sum(y == 0) == 0) return(NA_real_)
  sum(pred == 0 & y == 0) / sum(y == 0)
}

n_bm <- length(y_bm)
boot_mat <- matrix(NA_real_, nrow = B, ncol = 6,
                   dimnames = list(NULL, c("AUROC", "AUPRC", "Brier",
                                           "ECE", "Sensitivity", "Specificity")))
for (b in seq_len(B)) {
  idx         <- sample(n_bm, n_bm, replace = TRUE)
  yb <- y_bm[idx]; pb <- p_bm[idx]
  boot_mat[b, "AUROC"]       <- auroc(yb, pb)
  boot_mat[b, "AUPRC"]       <- auprc(yb, pb)
  boot_mat[b, "Brier"]       <- brier(yb, pb)
  boot_mat[b, "ECE"]         <- ece(yb, pb)
  boot_mat[b, "Sensitivity"] <- sens_at_half(yb, pb)
  boot_mat[b, "Specificity"] <- spec_at_half(yb, pb)
}

point_ests <- c(
  AUROC       = auroc(y_bm, p_bm),
  AUPRC       = auprc(y_bm, p_bm),
  Brier       = brier(y_bm, p_bm),
  ECE         = ece(y_bm, p_bm),
  Sensitivity = sens_at_half(y_bm, p_bm),
  Specificity = spec_at_half(y_bm, p_bm)
)

boot_ci <- do.call(rbind, lapply(colnames(boot_mat), function(m) {
  vals <- boot_mat[, m]
  ci   <- quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    metric    = m,
    estimate  = point_ests[[m]],
    ci_lower  = ci[1],
    ci_upper  = ci[2],
    B         = B
  )
}))

utils::write.csv(boot_ci, "outputs/tables/bootstrapped_metrics_best_model.csv",
                 row.names = FALSE)
message("[08b] Bootstrapped CIs saved.")

## =============================================================================
## 6. NRI AND IDI — BEST MODEL vs LOGISTIC REGRESSION BASELINE
## =============================================================================

baseline_name <- "glmnet"

p_base <- preds_all[[baseline_name]]
ok_nri <- !is.na(p_bm) & !is.na(p_base) & !is.na(y_all)

y_nri   <- y_all[ok_nri]
p_new   <- p_bm[ok_nri]
p_old   <- p_base[ok_nri]

## point estimate NRI (continuous, no categories)
nri_point <- function(y, p_new, p_old) {
  events     <- y == 1
  non_events <- y == 0
  ## NRI for events: proportion correctly moved up (new > old) minus down
  up_ev   <- mean(p_new[events]     > p_old[events])
  down_ev <- mean(p_new[events]     < p_old[events])
  up_nev  <- mean(p_new[non_events] > p_old[non_events])
  down_nev<- mean(p_new[non_events] < p_old[non_events])
  nri_ev  <- up_ev - down_ev
  nri_nev <- down_nev - up_nev   # correct for non-events: move DOWN is good
  list(NRI_events     = nri_ev,
       NRI_non_events = nri_nev,
       NRI_total      = nri_ev + nri_nev)
}

## IDI
idi_point <- function(y, p_new, p_old) {
  events     <- y == 1
  non_events <- y == 0
  is1 <- mean(p_new[events])     - mean(p_old[events])
  is0 <- mean(p_new[non_events]) - mean(p_old[non_events])
  list(IDI = is1 - is0, IS1 = is1, IS0 = is0)
}

nri_est <- nri_point(y_nri, p_new, p_old)
idi_est <- idi_point(y_nri, p_new, p_old)

## bootstrap CIs
set.seed(20260427)
n_nri     <- length(y_nri)
nri_boot  <- matrix(NA_real_, B, 4,
                    dimnames = list(NULL, c("NRI_events","NRI_non_events",
                                            "NRI_total","IDI")))
for (b in seq_len(B)) {
  idx  <- sample(n_nri, n_nri, replace = TRUE)
  yb   <- y_nri[idx]; pnb <- p_new[idx]; pob <- p_old[idx]
  if (length(unique(yb)) < 2) next
  nr   <- nri_point(yb, pnb, pob)
  id   <- idi_point(yb, pnb, pob)
  nri_boot[b, ] <- c(nr$NRI_events, nr$NRI_non_events, nr$NRI_total, id$IDI)
}

nri_idi_df <- data.frame(
  metric     = c("NRI_events", "NRI_non_events", "NRI_total", "IDI"),
  estimate   = c(nri_est$NRI_events, nri_est$NRI_non_events,
                 nri_est$NRI_total, idi_est$IDI),
  ci_lower   = apply(nri_boot, 2, quantile, 0.025, na.rm = TRUE),
  ci_upper   = apply(nri_boot, 2, quantile, 0.975, na.rm = TRUE),
  model_new  = best_name,
  model_base = baseline_name,
  B          = B
)

utils::write.csv(nri_idi_df, "outputs/tables/nri_idi_vs_baseline.csv",
                 row.names = FALSE)
message("[08b] NRI/IDI saved.")

## =============================================================================
## 7. PERFORMANCE BY CONFIDENCE TIER
## =============================================================================

ok_tier <- !is.na(p_bm) & !is.na(y_bm)
y_tier  <- y_bm[ok_tier]
p_tier  <- p_bm[ok_tier]
n_test  <- length(y_tier)

tier <- ifelse(p_tier > 0.80, "High-confidence P/LP",
               ifelse(p_tier < 0.20, "High-confidence B/LB", "Uncertain/Abstain"))

tier_rows <- lapply(c("High-confidence P/LP", "High-confidence B/LB",
                       "Uncertain/Abstain"), function(t) {
  idx <- which(tier == t)
  n_t <- length(idx)
  if (n_t == 0) return(data.frame(
    tier            = t, n = 0, prevalence = NA_real_,
    proportion      = 0, AUROC_within = NA_real_,
    PPV_boundary    = NA_real_, NPV_boundary = NA_real_
  ))

  y_t <- y_tier[idx]; p_t <- p_tier[idx]
  prev_t <- mean(y_t)

  ## AUROC within tier (only meaningful for high-confidence combined)
  auroc_t <- if (length(unique(y_t)) > 1) auroc(y_t, p_t) else NA_real_

  ## PPV and NPV at tier boundaries
  ppv <- npv <- NA_real_
  if (t == "High-confidence P/LP") {
    ## boundary = 0.80; all in tier are predicted positive
    ppv <- mean(y_t)            # all predicted P/LP; PPV = obs prevalence
    ## NPV: among those below boundary (≤ 0.80) — computed across full set
    below <- which(p_tier <= 0.80)
    if (length(below) > 0) npv <- mean(y_tier[below] == 0)
  } else if (t == "High-confidence B/LB") {
    ## boundary = 0.20; all in tier predicted negative
    npv <- mean(y_t == 0)
    above <- which(p_tier >= 0.20)
    if (length(above) > 0) ppv <- mean(y_tier[above])
  }

  data.frame(
    tier            = t,
    n               = n_t,
    prevalence      = prev_t,
    proportion      = n_t / n_test,
    AUROC_within    = auroc_t,
    PPV_boundary    = ppv,
    NPV_boundary    = npv
  )
})

## combined high-confidence tier AUROC
hc_idx  <- which(tier %in% c("High-confidence P/LP", "High-confidence B/LB"))
auroc_hc <- if (length(unique(y_tier[hc_idx])) > 1)
  auroc(y_tier[hc_idx], p_tier[hc_idx]) else NA_real_

tier_df <- do.call(rbind, tier_rows)
tier_df <- rbind(tier_df,
                 data.frame(
                   tier         = "High-confidence (combined)",
                   n            = length(hc_idx),
                   prevalence   = mean(y_tier[hc_idx]),
                   proportion   = length(hc_idx) / n_test,
                   AUROC_within = auroc_hc,
                   PPV_boundary = NA_real_,
                   NPV_boundary = NA_real_
                 ))

utils::write.csv(tier_df, "outputs/tables/performance_by_confidence_tier.csv",
                 row.names = FALSE)
message("[08b] Performance by confidence tier saved.")

## =============================================================================
## WRAP-UP
## =============================================================================
message(
  "\n[08b] All outputs written:\n",
  "  Tables:\n",
  "    outputs/tables/calibration_metrics_all_learners.csv\n",
  "    outputs/tables/dca_results.csv\n",
  "    outputs/tables/delong_pairwise_tests.csv\n",
  "    outputs/tables/bootstrapped_metrics_best_model.csv\n",
  "    outputs/tables/nri_idi_vs_baseline.csv\n",
  "    outputs/tables/performance_by_confidence_tier.csv\n",
  "  Figures:\n",
  "    outputs/figures/reliability_diagram_best_model.png\n",
  "    outputs/figures/reliability_diagrams_all_learners.png\n",
  "    outputs/figures/decision_curve_analysis.png\n"
)
