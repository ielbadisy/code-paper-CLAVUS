# ==============================================================================
# 12_pdp_ice_figure.R  —  PDP + ICE figure using correct top-5 VIP features
#
# The saved interpretability object used hardcoded external annotation scores
# as "top features"; this script recomputes PDP and ICE for the actual top-5
# features by XGBoost gain importance and produces Figure 9.
# ==============================================================================

source("R/00_setup.R")
library(funcml)
library(ggplot2)
library(patchwork)

# ── 0. Okabe-Ito palette ────────────────────────────────────────────────────
OI <- c(black="#000000", orange="#E69F00", sky_blue="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", pink="#CC79A7",
        grey="#999999")

FEAT_LABELS <- c(
  days_since_first_seen        = "Days since first ClinVar observation",
  gene_pathogenic_fraction_at_t0 = "Gene pathogenic fraction at baseline",
  gene_benign_count_at_t0      = "Gene benign variant count at baseline",
  baseline_number_submitters   = "Number of ClinVar submitters",
  baseline_stars_numeric       = "ClinVar review star rating"
)

# ── 1. Load data and model ───────────────────────────────────────────────────
dat  <- readRDS("data/curated/vus_resolution_dataset.rds")
best <- readRDS("outputs/models/best_model.rds")
imp_bundle <- readRDS("outputs/models/mimar_rf_imputer.rds")

# ── 2. Prepare funcml data (same as 09_interpretability.R) ──────────────────
clean_funcml_frame <- function(df, predictors) {
  out <- df[, c("label", predictors), drop = FALSE]
  for (nm in predictors) {
    if (is.logical(out[[nm]])) out[[nm]] <- as.integer(out[[nm]])
    if (is.character(out[[nm]])) out[[nm]][is.na(out[[nm]])] <- "Missing"
    if (is.numeric(out[[nm]])) {
      med <- stats::median(out[[nm]], na.rm = TRUE)
      if (is.na(med)) med <- 0
      out[[nm]][is.na(out[[nm]])] <- med
    }
  }
  if (length(best$funcml_levels)) {
    for (nm in names(best$funcml_levels)) {
      vals <- as.character(out[[nm]])
      vals[!vals %in% best$funcml_levels[[nm]]] <- "OTHER"
      out[[nm]] <- factor(vals, levels = best$funcml_levels[[nm]])
    }
  }
  out$label <- factor(ifelse(out$label == 1, "P_LP", "B_LB"), levels = c("B_LB", "P_LP"))
  out
}

funcml_data_raw <- clean_funcml_frame(dat[dat$split_temporal == "test", ], best$predictors)
feat09 <- funcml_data_raw[, best$predictors, drop = FALSE]

for (col in imp_bundle$cols_to_imp) {
  if (is.null(imp_bundle$models[[col]])) next
  miss <- is.na(feat09[[col]]); if (!any(miss)) next
  others <- setdiff(imp_bundle$numeric_preds, col)
  others <- others[others %in% names(feat09)]
  x_new <- feat09[miss, others, drop = FALSE]
  for (oc in others) {
    med <- median(feat09[[oc]], na.rm = TRUE)
    x_new[[oc]][is.na(x_new[[oc]])] <- if (is.na(med)) 0 else med
  }
  tryCatch(
    feat09[[col]][miss] <- stats::predict(imp_bundle$models[[col]], newdata = x_new),
    error = function(e) NULL
  )
}
for (col in best$predictors) {
  if (is.numeric(feat09[[col]]) && any(is.na(feat09[[col]]))) {
    med <- median(feat09[[col]], na.rm = TRUE)
    feat09[[col]][is.na(feat09[[col]])] <- if (is.na(med)) 0 else med
  }
}
funcml_data <- cbind(label = funcml_data_raw$label, feat09)
set.seed(20260427)
if (nrow(funcml_data) > 1000) funcml_data <- funcml_data[sample(seq_len(nrow(funcml_data)), 1000), ]

# ── 3. Top 5 continuous/ordinal features from actual VIP ────────────────────
top5 <- c("days_since_first_seen", "gene_pathogenic_fraction_at_t0",
           "gene_benign_count_at_t0", "baseline_number_submitters",
           "baseline_stars_numeric")
top5 <- intersect(top5, best$predictors)

# ── 4. Compute PDP via funcml::interpret ────────────────────────────────────
message("[12] Computing PDP for: ", paste(top5, collapse = ", "))
pdp_obj <- tryCatch(
  funcml::interpret(fit = best$model, data = funcml_data,
                    method = "pdp", features = top5,
                    type = "prob", class_level = "P_LP"),
  error = function(e) { message("PDP error: ", e$message); NULL }
)

# ── 5. Compute ICE via funcml::interpret (top 3 continuous) ─────────────────
ice_feats <- c("days_since_first_seen", "gene_pathogenic_fraction_at_t0",
               "baseline_number_submitters")
ice_feats <- intersect(ice_feats, best$predictors)
message("[12] Computing ICE for: ", paste(ice_feats, collapse = ", "))
ice_obj <- tryCatch(
  funcml::interpret(fit = best$model, data = funcml_data,
                    method = "ice", features = ice_feats,
                    type = "prob", class_level = "P_LP"),
  error = function(e) { message("ICE error: ", e$message); NULL }
)

# ── 6. Extract curve data ────────────────────────────────────────────────────
pdp_df <- if (!is.null(pdp_obj) && "result" %in% names(pdp_obj))
  pdp_obj$result$curves else NULL
ice_df <- if (!is.null(ice_obj) && "result" %in% names(ice_obj))
  ice_obj$result$curves else NULL

# ── 7. Build plots ───────────────────────────────────────────────────────────
theme_pub <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 10, face = "bold"))

## Panel A: PDP
if (!is.null(pdp_df)) {
  pdp_df$feat_label <- FEAT_LABELS[pdp_df$feature]
  pdp_df$feat_label[is.na(pdp_df$feat_label)] <- pdp_df$feature[is.na(pdp_df$feat_label)]
  pdp_df$feat_label <- factor(pdp_df$feat_label, levels = FEAT_LABELS)

  pA <- ggplot(pdp_df, aes(x = value, y = yhat)) +
    geom_line(colour = OI[["blue"]], linewidth = 0.9) +
    geom_point(colour = OI[["blue"]], size = 1.8) +
    facet_wrap(~ feat_label, scales = "free_x", ncol = 3) +
    labs(title = "A  Partial dependence plots (top 5 features)",
         x = "Feature value", y = "Predicted P/LP probability (PDP)") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pub
} else {
  pA <- ggplot() + annotate("text", x=0, y=0, label="PDP failed") + theme_void()
}

## Panel B: ICE
if (!is.null(ice_df)) {
  ice_df$feat_label <- FEAT_LABELS[ice_df$feature]
  ice_df$feat_label[is.na(ice_df$feat_label)] <- ice_df$feature[is.na(ice_df$feat_label)]
  # compute per-feature PDP (mean) for overlay
  ice_mean <- aggregate(yhat ~ feature + value + feat_label, data = ice_df, FUN = mean)

  pB <- ggplot(ice_df, aes(x = value, y = yhat, group = id)) +
    geom_line(colour = "grey70", alpha = 0.25, linewidth = 0.3) +
    geom_line(data = ice_mean, aes(group = NULL),
              colour = OI[["vermillion"]], linewidth = 1.2) +
    facet_wrap(~ feat_label, scales = "free_x", ncol = 3) +
    labs(title = "B  Individual conditional expectation curves (ICE; bold = PDP average)",
         x = "Feature value", y = "Predicted P/LP probability") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pub
} else {
  pB <- ggplot() + annotate("text", x=0, y=0, label="ICE failed") + theme_void()
}

# ── 8. Combine and save ──────────────────────────────────────────────────────
fig9 <- pA / pB + plot_layout(heights = c(1, 1))

ragg::agg_png("outputs/figures/figure9_pdp_ice_curves.png",
              width = 183, height = 200, units = "mm", res = 400)
print(fig9)
dev.off()

cairo_pdf("outputs/figures/figure9_pdp_ice_curves.pdf",
          width = 183/25.4, height = 200/25.4)
print(fig9)
dev.off()

message("[12] Figure 9 saved: outputs/figures/figure9_pdp_ice_curves.{png,pdf}")
