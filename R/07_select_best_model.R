source("R/00_setup.R")

performance <- read.csv("outputs/tables/model_performance.csv")
performance$slope_distance <- abs(performance$calibration_slope - 1)
simplicity_order <- c(
  "logistic_regression",
  "penalized_logistic_regression",
  "decision_tree",
  "mars_earth",
  "c50_rules",
  "radial_svm",
  "neural_net",
  "random_forest",
  "gbm_boosting",
  "gradient_boosting",
  "lightgbm_boosting",
  "stacked_ensemble"
)
simplicity <- seq_along(simplicity_order)
names(simplicity) <- simplicity_order
performance$simplicity <- ifelse(performance$learner %in% names(simplicity), simplicity[performance$learner], length(simplicity) + 1)
ord <- order(-performance$AUPRC, performance$Brier_score, performance$slope_distance, performance$simplicity)
best_name <- performance$learner[ord[1]]

bundle <- readRDS("outputs/models/all_models.rds")
best_model <- list(
  name = best_name,
  model = bundle$models[[best_name]],
  all_models = bundle$models,
  calibrators = bundle$calibrators,
  imputer = bundle$imputer,
  predictors = bundle$predictors,
  funcml_levels = bundle$funcml_levels
)

saveRDS(best_model, "outputs/models/best_model.rds")
writeLines(best_name, "outputs/models/best_model_name.txt")
utils::write.csv(performance[ord, ], "outputs/tables/model_selection_ranking.csv", row.names = FALSE)
log_message("outputs/logs/pipeline_warnings.txt", sprintf("Selected best model by temporal-test AUPRC hierarchy: %s.", best_name))
