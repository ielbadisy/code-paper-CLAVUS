source("R/00_setup.R")
load_or_stop("funcml")
library(mimar)
load_or_stop("ggplot2")

dat <- readRDS("data/curated/vus_resolution_dataset.rds")
best <- readRDS("outputs/models/best_model.rds")
pred <- read.csv("outputs/predictions/best_model_temporal_test_predictions.csv")

perf <- read.csv("outputs/tables/model_performance.csv")
model_name <- best$name

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
  if (exists("best") && length(best$funcml_levels)) {
    for (nm in names(best$funcml_levels)) {
      vals <- as.character(out[[nm]])
      vals[!vals %in% best$funcml_levels[[nm]]] <- "OTHER"
      out[[nm]] <- factor(vals, levels = best$funcml_levels[[nm]])
    }
  }
  out$label <- factor(ifelse(out$label == 1, "P_LP", "B_LB"), levels = c("B_LB", "P_LP"))
  out
}

importance_from_model <- function(funcml_permute = NULL) {
  if (!is.null(funcml_permute) && inherits(funcml_permute, "funcml_permute")) {
    candidate <- funcml_permute
    if ("result" %in% names(candidate) && "scores" %in% names(candidate$result) && is.data.frame(candidate$result$scores)) {
      df <- candidate$result$scores
    } else if ("importance" %in% names(candidate) && is.data.frame(candidate$importance)) {
      df <- candidate$importance
    } else if ("results" %in% names(candidate) && is.data.frame(candidate$results)) {
      df <- candidate$results
    } else {
      df <- NULL
    }
    if (!is.null(df)) {
      v <- intersect(c("feature", "variable", "term"), names(df))[1]
      s <- intersect(c("importance", "dropout_loss", "value", "mean"), names(df))[1]
      if (!is.na(v) && !is.na(s)) return(data.frame(variable = df[[v]], importance = abs(df[[s]])))
    }
  }
  if (model_name == "logistic_regression") {
    co <- coef(best$model$state$state)
    data.frame(variable = names(co)[-1], importance = abs(as.numeric(co[-1])))
  } else if (model_name == "super_learner_stacked_ensemble") {
    co <- coef(best$model)
    data.frame(variable = names(co)[-1], importance = abs(as.numeric(co[-1])))
  } else if (model_name == "penalized_logistic_regression") {
    co <- as.matrix(coef(best$model$state$state))
    data.frame(variable = rownames(co)[-1], importance = abs(as.numeric(co[-1, 1])))
  } else if (model_name == "random_forest") {
    data.frame(variable = best$predictors, importance = NA_real_)
  } else {
    imp <- xgboost::xgb.importance(model = best$model$state$state)
    data.frame(variable = imp$Feature, importance = imp$Gain)
  }
}

funcml_fit_for_interpret <- best$model
funcml_permute <- NULL
funcml_calibration <- NULL
interpretability_results <- list()
interpretability_status <- list()
if (inherits(funcml_fit_for_interpret, "funcml_fit")) {
  funcml_data_raw <- clean_funcml_frame(dat[dat$split_temporal == "test", ], best$predictors)
  imp_bundle <- readRDS("outputs/models/mimar_rf_imputer.rds")
  feat09 <- funcml_data_raw[, best$predictors, drop = FALSE]
  for (col in imp_bundle$cols_to_imp) {
    if (is.null(imp_bundle$models[[col]])) next
    miss <- is.na(feat09[[col]]); if (!any(miss)) next
    others <- setdiff(imp_bundle$numeric_preds, col); others <- others[others %in% names(feat09)]
    x_new <- feat09[miss, others, drop=FALSE]
    for (oc in others) { med <- median(feat09[[oc]], na.rm=TRUE); x_new[[oc]][is.na(x_new[[oc]])] <- if(is.na(med)) 0 else med }
    tryCatch(feat09[[col]][miss] <- stats::predict(imp_bundle$models[[col]], newdata=x_new), error=function(e) NULL)
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

  top_features <- c("gnomad_af", "cadd_score", "revel_score", "alphamissense_score",
                    "baseline_stars_numeric", "gene_pathogenic_fraction_at_t0",
                    "gene_benign_fraction_at_t0", "baseline_number_submitters")
  top_features <- intersect(top_features, best$predictors)
  if (!length(top_features)) top_features <- head(best$predictors, 5)

  method_specs <- list(
    vip = list(args = list(method = "vip", type = "prob", class_level = "P_LP")),
    permute = list(args = list(method = "permute", type = "prob", metric = "brier", class_level = "P_LP", seed = 20260427)),
    pdp = list(args = list(method = "pdp", features = head(top_features, 5), type = "prob", class_level = "P_LP")),
    ice = list(args = list(method = "ice", features = head(top_features, 3), type = "prob", class_level = "P_LP")),
    ale = list(args = list(method = "ale", features = head(top_features, 5), type = "prob", class_level = "P_LP")),
    local = list(args = list(method = "local", newdata = head(funcml_data, 3), type = "prob", class_level = "P_LP")),
    lime = list(args = list(method = "lime", newdata = head(funcml_data, 3), type = "prob", class_level = "P_LP")),
    shap = list(args = list(method = "shap", newdata = head(funcml_data, 3), nsim = 25, type = "prob", class_level = "P_LP")),
    local_model = list(args = list(method = "local_model", newdata = head(funcml_data, 3), type = "prob", class_level = "P_LP")),
    interaction = list(args = list(method = "interaction", features = head(top_features, 5), type = "prob", class_level = "P_LP")),
    surrogate = list(args = list(method = "surrogate", type = "prob", class_level = "P_LP")),
    profile = list(args = list(method = "profile", features = head(top_features, 5), type = "prob", class_level = "P_LP")),
    ceteris_paribus = list(args = list(method = "ceteris_paribus", newdata = head(funcml_data, 3), features = head(top_features, 5), type = "prob", class_level = "P_LP")),
    calibration = list(args = list(method = "calibration", type = "prob", class_level = "P_LP", bins = 8, strategy = "quantile"))
  )

  for (nm in names(method_specs)) {
    res <- tryCatch(
      do.call(funcml::interpret, c(list(fit = funcml_fit_for_interpret, data = funcml_data), method_specs[[nm]]$args)),
      error = function(e) e
    )
    ok <- !inherits(res, "error")
    interpretability_status[[nm]] <- data.frame(method = nm, succeeded = ok, note = if (ok) "ok" else conditionMessage(res))
    interpretability_results[[nm]] <- res
    if (ok) {
      png(file.path("outputs/figures", paste0("funcml_interpret_", nm, "_best_model.png")), width = 1100, height = 800)
      plot_res <- try(print(plot(res)), silent = TRUE)
      dev.off()
      if (inherits(plot_res, "try-error")) interpretability_status[[nm]]$note <- paste("computed; plot failed", as.character(plot_res))
    }
  }
  funcml_permute <- interpretability_results$permute
  funcml_calibration <- interpretability_results$calibration
  saveRDS(interpretability_results, "outputs/models/funcml_interpretability_best_model.rds")
  utils::write.csv(do.call(rbind, interpretability_status), "outputs/tables/funcml_interpretability_status.csv", row.names = FALSE)
  if (inherits(funcml_permute, "funcml_permute")) {
    file.copy("outputs/figures/funcml_interpret_permute_best_model.png", "outputs/figures/funcml_permutation_importance_best_model.png", overwrite = TRUE)
  }
  if (inherits(funcml_calibration, "funcml_calibration")) {
    file.copy("outputs/figures/funcml_interpret_calibration_best_model.png", "outputs/figures/funcml_calibration_best_model.png", overwrite = TRUE)
  }
}

vi <- importance_from_model(funcml_permute)
vi <- vi[!is.na(vi$importance), ]
vi <- vi[order(vi$importance, decreasing = TRUE), ]
utils::write.csv(vi, "outputs/tables/variable_importance_best_model.csv", row.names = FALSE)

png("outputs/figures/variable_importance_best_model.png", width = 1100, height = 800)
print(ggplot2::ggplot(head(vi, 20), ggplot2::aes(x = reorder(variable, importance), y = importance)) +
        ggplot2::geom_col(fill = "#2f6f73") + ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Importance", title = "Variable importance for selected best model") +
        ggplot2::theme_minimal(base_size = 13))
dev.off()

cal <- aggregate(cbind(label, best_model_score) ~ cut(best_model_score, seq(0, 1, .1), include.lowest = TRUE), pred, mean)
names(cal)[1] <- "bin"
png("outputs/figures/calibration_curve_best_model.png", width = 900, height = 700)
print(ggplot2::ggplot(cal, ggplot2::aes(x = best_model_score, y = label)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray40") +
        ggplot2::geom_point(size = 3, color = "#8b3a3a") + ggplot2::geom_line(color = "#8b3a3a") +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "Mean predicted probability", y = "Observed P/LP resolution", title = "Calibration curve") +
        ggplot2::theme_minimal(base_size = 13))
dev.off()

png("outputs/figures/prediction_distribution_by_class.png", width = 950, height = 650)
print(ggplot2::ggplot(pred, ggplot2::aes(x = best_model_score, fill = final_class_group)) +
        ggplot2::geom_density(alpha = .45) +
        ggplot2::labs(x = "Predicted probability of later P/LP resolution", y = "Density", fill = "Final class") +
        ggplot2::theme_minimal(base_size = 13))
dev.off()

png("outputs/figures/risk_group_distribution.png", width = 900, height = 650)
print(ggplot2::ggplot(pred, ggplot2::aes(x = risk_group, fill = final_class_group)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::labs(x = NULL, y = "Proportion", fill = "Final class") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)))
dev.off()
