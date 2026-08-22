source("R/00_setup.R")
load_or_stop("funcml")
library(mimar)

dat <- readRDS("data/curated/vus_resolution_dataset.rds")
best <- readRDS("outputs/models/best_model.rds")
imp_bundle <- readRDS("outputs/models/mimar_rf_imputer.rds")

prepare_rows <- function(rows) {
  predictors <- best$predictors
  out <- rows[, c("label", predictors), drop = FALSE]
  for (nm in predictors) {
    if (is.logical(out[[nm]])) out[[nm]] <- as.integer(out[[nm]])
    if (is.character(out[[nm]])) out[[nm]][is.na(out[[nm]])] <- "Missing"
  }
  for (nm in names(best$funcml_levels)) {
    vals <- as.character(out[[nm]])
    vals[is.na(vals)] <- "Missing"
    vals[!vals %in% best$funcml_levels[[nm]]] <- "OTHER"
    out[[nm]] <- factor(vals, levels = best$funcml_levels[[nm]])
  }
  # Apply pre-fitted mimar RF imputer (per-variable)
  feat <- out[, predictors, drop = FALSE]
  for (col in imp_bundle$cols_to_imp) {
    if (is.null(imp_bundle$models[[col]])) next
    miss <- is.na(feat[[col]])
    if (!any(miss)) next
    others <- setdiff(imp_bundle$numeric_preds, col)
    others <- others[others %in% names(feat)]
    x_new  <- feat[miss, others, drop = FALSE]
    for (oc in others) { med <- median(feat[[oc]], na.rm=TRUE); x_new[[oc]][is.na(x_new[[oc]])] <- if(is.na(med)) 0 else med }
    tryCatch(feat[[col]][miss] <- stats::predict(imp_bundle$models[[col]], newdata=x_new),
             error=function(e) NULL)
  }
  # Final safety pass: median-fill any remaining NAs in numeric columns
  for (col in predictors) {
    if (is.numeric(feat[[col]]) && any(is.na(feat[[col]]))) {
      med <- median(feat[[col]], na.rm = TRUE)
      feat[[col]][is.na(feat[[col]])] <- if (is.na(med)) 0 else med
    }
  }
  cbind(label = factor(ifelse(out$label == 1, "P_LP", "B_LB"), levels = c("B_LB", "P_LP")), feat)
}

predict_base <- function(model_name, rows) {
  model_rows <- prepare_rows(rows)
  m <- best$all_models[[model_name]]
  safe_prob(stats::predict(m, newdata = model_rows, type = "prob")[, "P_LP"])
}

predict_named <- function(model_name, rows) {
  p <- predict_base(model_name, rows)
  safe_prob(stats::predict(best$calibrators[[model_name]], newdata = data.frame(p = p), type = "response"))
}

metric_row <- function(y, p, validation_scheme) {
  pred <- as.integer(p >= 0.5)
  data.frame(
    validation_scheme = validation_scheme,
    n = length(y),
    AUROC = if (length(unique(y)) < 2) NA_real_ else {
      r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    },
    AUPRC = {
      ord <- order(p, decreasing = TRUE); yy <- y[ord]; tp <- cumsum(yy == 1); fp <- cumsum(yy == 0)
      sum(diff(c(0, tp / max(sum(yy == 1), 1))) * (tp / pmax(tp + fp, 1)))
    },
    Brier_score = mean((y - p)^2),
    sensitivity = sum(pred == 1 & y == 1) / max(sum(y == 1), 1),
    specificity = sum(pred == 0 & y == 0) / max(sum(y == 0), 1)
  )
}

temporal <- dat[dat$split_temporal == "test", ]
gene <- dat[dat$split_gene_heldout == "test", ]
temporal$best_model_score <- predict_named(best$name, temporal)
gene$best_model_score <- predict_named(best$name, gene)

validation <- rbind(
  metric_row(temporal$label, temporal$best_model_score, "temporal_test"),
  metric_row(gene$label, gene$best_model_score, "gene_held_out")
)

temporal$risk_group <- cut(temporal$best_model_score, c(-Inf, .2, .8, Inf),
                           labels = c("likely B/LB resolution", "remain uncertain / abstain", "likely P/LP resolution"))
utils::write.csv(validation, "outputs/tables/best_model_validation.csv", row.names = FALSE)
utils::write.csv(temporal[, c("variant_id", "label", "final_class_group", "best_model_score", "risk_group")],
                 "outputs/predictions/best_model_temporal_test_predictions.csv", row.names = FALSE)

subgroup_metric <- function(rows, group_name, group_value) {
  if (nrow(rows) < 50 || length(unique(rows$label)) < 2) return(NULL)
  cbind(group = group_name, level = as.character(group_value),
        metric_row(rows$label, rows$best_model_score, paste0("subgroup_", group_name)))
}
subgroups <- list()
for (v in c("final_release_month", "baseline_review_status", "variant_type")) {
  if (v %in% names(temporal)) {
    pieces <- lapply(split(temporal, temporal[[v]]), function(z) subgroup_metric(z, v, z[[v]][1]))
    subgroups <- c(subgroups, pieces[!vapply(pieces, is.null, logical(1))])
  }
}
if (length(subgroups)) {
  utils::write.csv(do.call(rbind, subgroups), "outputs/tables/sensitivity_subgroup_validation.csv", row.names = FALSE)
}
