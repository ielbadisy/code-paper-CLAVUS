source("R/00_setup.R")
library(funcml)
library(mimar)

dat     <- readRDS("data/curated/vus_resolution_dataset.rds")
recipe  <- readRDS("data/curated/feature_recipe.rds")
predictors <- recipe$predictor_names

# ── Splits ────────────────────────────────────────────────────────────────────
train_full <- dat[dat$split_temporal == "train", ]
test       <- dat[dat$split_temporal == "test",  ]
train_full <- train_full[order(train_full$final_date, train_full$variant_id), ]
cal_start  <- max(2, floor(0.85 * nrow(train_full)))
calibration <- train_full[seq(cal_start, nrow(train_full)), ]
train       <- train_full[seq_len(cal_start - 1), ]

# ── Drop constant predictors ─────────────────────────────────────────────────
predictors <- predictors[vapply(predictors, function(nm)
  length(unique(train[[nm]][!is.na(train[[nm]])])) >= 2, logical(1))]

# ── Prep raw frames (factor levels, NA→"Missing" for chars) ──────────────────
prepare_raw <- function(df, predictors) {
  out <- df[, c("label", predictors), drop = FALSE]
  for (nm in predictors) {
    if (is.logical(out[[nm]])) out[[nm]] <- as.integer(out[[nm]])
    if (is.character(out[[nm]])) out[[nm]][is.na(out[[nm]])] <- "Missing"
  }
  out
}

raw_train <- prepare_raw(train,       predictors)
raw_cal   <- prepare_raw(calibration, predictors)
raw_test  <- prepare_raw(test,        predictors)

# ── Factor level alignment ────────────────────────────────────────────────────
funcml_levels <- list()
for (nm in predictors) {
  if (is.character(raw_train[[nm]]) || is.factor(raw_train[[nm]])) {
    lv <- unique(c(as.character(raw_train[[nm]]), "OTHER", "Missing"))
    raw_train[[nm]] <- factor(as.character(raw_train[[nm]]), levels = lv)
    align <- function(x, lv) { x <- as.character(x); x[!x %in% lv] <- "OTHER"; factor(x, levels = lv) }
    raw_cal[[nm]]  <- align(raw_cal[[nm]],  lv)
    raw_test[[nm]] <- align(raw_test[[nm]], lv)
    funcml_levels[[nm]] <- lv
  }
}

# ── mimar RF imputation: fit on training, apply to cal/test ──────────────────
numeric_preds <- predictors[vapply(raw_train[, predictors, drop = FALSE], is.numeric, logical(1))]
cols_to_imp   <- numeric_preds[vapply(raw_train[, numeric_preds, drop = FALSE],
                                      function(x) any(is.na(x)), logical(1))]

message("[06] Fitting mimar RF imputers for ", length(cols_to_imp), " variables ...")

imp_models <- vector("list", length(cols_to_imp))
names(imp_models) <- cols_to_imp

for (col in cols_to_imp) {
  obs  <- !is.na(raw_train[[col]])
  if (sum(obs) < 10) next
  others <- setdiff(numeric_preds, col)
  others <- others[vapply(raw_train[obs, others, drop = FALSE],
                          function(x) !all(is.na(x)), logical(1))]
  if (length(others) == 0) next
  x_obs <- raw_train[obs, others, drop = FALSE]
  for (oc in others) {                              # median-fill predictor block NAs
    med <- stats::median(x_obs[[oc]], na.rm = TRUE)
    x_obs[[oc]][is.na(x_obs[[oc]])] <- if (is.na(med)) 0 else med
  }
  tryCatch({
    obj <- mimar::imputer("rf")
    imp_models[[col]] <- mimar::fit(obj, x = x_obs, y = raw_train[[col]][obs],
                                    seed = 20260427)
  }, error = function(e) message("  mimar fit skipped for ", col, ": ", e$message))
}

apply_imputations <- function(df, imp_models, numeric_preds) {
  for (col in names(imp_models)) {
    if (is.null(imp_models[[col]])) next
    miss <- is.na(df[[col]])
    if (!any(miss)) next
    others <- setdiff(numeric_preds, col)
    others <- others[others %in% names(df)]
    x_new <- df[miss, others, drop = FALSE]
    for (oc in others) {
      med <- stats::median(df[[oc]], na.rm = TRUE)
      x_new[[oc]][is.na(x_new[[oc]])] <- if (is.na(med)) 0 else med
    }
    tryCatch({
      df[[col]][miss] <- stats::predict(imp_models[[col]], newdata = x_new)
    }, error = function(e) {
      med <- stats::median(df[[col]], na.rm = TRUE)
      df[[col]][miss] <- if (is.na(med)) 0 else med
    })
  }
  df
}

funcml_train <- apply_imputations(raw_train, imp_models, numeric_preds)
funcml_cal   <- apply_imputations(raw_cal,   imp_models, numeric_preds)
funcml_test  <- apply_imputations(raw_test,  imp_models, numeric_preds)

# Outcome as factor
to_factor <- function(x) factor(ifelse(x == 1, "P_LP", "B_LB"), levels = c("B_LB", "P_LP"))
funcml_train$label <- to_factor(funcml_train$label)
funcml_cal$label   <- to_factor(funcml_cal$label)
funcml_test$label  <- to_factor(funcml_test$label)

y_train <- train$label
y_test  <- test$label
funcml_formula <- stats::as.formula(paste("label ~", paste(predictors, collapse = " + ")))

# ── Time-aware rolling CV for hyperparameter tuning ───────────────────────────
# funcml_train is already sorted chronologically; use row index as time ordering.
# Expanding window: 60% initial window, 20% assessment slices ≈ 2 folds.
n_train_rows   <- nrow(funcml_train)
time_resampling <- funcml::cv(
  method     = "time",
  time       = seq_len(n_train_rows),
  initial    = floor(n_train_rows * 0.60),
  assess     = floor(n_train_rows * 0.20),
  cumulative = TRUE,
  seed       = 20260427
)

# ── Display-name mapping (no underscores) ─────────────────────────────────────
DISPLAY <- c(
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

# ── Learner registry — all tuned ──────────────────────────────────────────────
learner_registry <- list(
  xgboost = list(model = "xgboost", mode = "tune",
    grid = expand.grid(nrounds   = c(100, 200),
                       max_depth = c(3, 5),
                       eta       = c(0.05, 0.10),
                       subsample = 0.8)),
  bart = list(model = "bart", mode = "tune",
    grid = expand.grid(ntrees = c(50, 100, 200))),
  ranger = list(model = "ranger", mode = "tune",
    grid = expand.grid(num.trees    = c(300, 500),
                       mtry         = c(4, 8),
                       min.node.size = c(5, 20))),
  gbm = list(model = "gbm", mode = "tune",
    grid = expand.grid(n.trees           = c(100, 200),
                       interaction.depth = c(2, 3),
                       shrinkage         = c(0.05, 0.10))),
  earth = list(model = "earth", mode = "tune",
    grid = expand.grid(degree = c(1, 2), nprune = c(10, 20))),
  glmnet = list(model = "glmnet", mode = "tune",
    grid = expand.grid(alpha  = c(0, 0.5, 1),
                       lambda = c(1e-4, 1e-3, 1e-2, 0.1))),
  svm = list(model = "e1071_svm", mode = "tune",
    grid = expand.grid(cost  = c(0.1, 1, 10),
                       gamma = c(0.01, 0.1))),
  nnet = list(model = "mlp", mode = "tune",
    grid  = expand.grid(hidden_units = c(32, 64, 128),
                        learn_rate   = c(0.01, 0.001)),
    extra = list(epochs = 300)),
  knn = list(model = "kknn", mode = "tune",
    grid = expand.grid(k = c(7, 15, 31)))
)

# ── Helper functions ──────────────────────────────────────────────────────────
predict_prob    <- function(fit, nd) safe_prob(stats::predict(fit, newdata = nd, type = "prob")[, "P_LP"])
fit_calibrator  <- function(y, p)   stats::glm(y ~ qlogis(safe_prob(p)), family = binomial())
apply_calibrator <- function(fit, p) safe_prob(stats::predict(fit, newdata = data.frame(p = safe_prob(p)), type = "response"))

auroc <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auprc <- function(y, p) {
  ord <- order(p, decreasing = TRUE); y <- y[ord]
  tp  <- cumsum(y == 1); fp <- cumsum(y == 0)
  pr  <- tp / pmax(tp + fp, 1); re <- tp / max(sum(y == 1), 1)
  sum(diff(c(0, re)) * pr)
}
compute_metrics <- function(y, p, nm) {
  p    <- safe_prob(p); pred <- as.integer(p >= 0.5)
  cal  <- tryCatch(stats::glm(y ~ qlogis(p), family = binomial()), error = function(e) NULL)
  # Use funcml built-in metrics; truth must be factor with positive class last
  truth_f <- factor(ifelse(y == 1, "P_LP", "B_LB"), levels = c("B_LB", "P_LP"))
  prob_mat <- cbind(B_LB = 1 - p, P_LP = p)
  data.frame(
    learner     = nm,
    display     = DISPLAY[nm],
    AUROC       = auroc(y, p),
    AUPRC       = auprc(y, p),
    Brier_score = funcml::brier(truth_f, prob_mat),
    log_loss    = funcml::logloss(truth_f, prob_mat),
    calibration_slope     = if (is.null(cal)) NA_real_ else unname(coef(cal)[2]),
    calibration_intercept = if (is.null(cal)) NA_real_ else unname(coef(cal)[1]),
    ECE         = tryCatch(funcml::ece(truth_f, prob_mat[, "P_LP"],
                                       positive = "P_LP"), error = function(e) NA_real_),
    sensitivity = sum(pred == 1 & y == 1) / max(sum(y == 1), 1),
    specificity = sum(pred == 0 & y == 0) / max(sum(y == 0), 1),
    PPV         = sum(pred == 1 & y == 1) / max(sum(pred == 1), 1),
    NPV         = sum(pred == 0 & y == 0) / max(sum(pred == 0), 1),
    F1_score    = {
      tp <- sum(pred == 1 & y == 1)
      2 * tp / max(2 * tp + sum(pred == 1 & y == 0) + sum(pred == 0 & y == 1), 1)
    },
    stringsAsFactors = FALSE
  )
}

# ── Train all learners ────────────────────────────────────────────────────────
models        <- list(); funcml_models <- list(); tuning_results <- list()
calibrators   <- list(); preds <- list(); raw_preds <- list()
calib_preds   <- list(); fit_status <- list(); best_params_list <- list()

for (nm in names(learner_registry)) {
  spec <- learner_registry[[nm]]
  message("[06] Tuning learner: ", DISPLAY[nm], " (", nm, ") ...")
  fitted <- tryCatch({
    tune_args <- c(
      list(funcml_train, funcml_formula,
           model      = spec$model,
           grid       = spec$grid,
           resampling = time_resampling,
           metric     = "brier",
           seed       = 20260427,
           verbose    = 0),
      if (!is.null(spec$extra)) spec$extra else list()
    )
    tuned <- do.call(funcml::tune, tune_args)
    tuning_results[[nm]] <- tuned
    best_params_list[[nm]] <- if (!is.null(tuned$best_params)) tuned$best_params else data.frame()
    tuned$fit_best
  }, error = function(e) {
    message("  Tuning failed: ", e$message, " — fitting with defaults.")
    tryCatch({
      fit_args <- c(
        list(funcml_formula, data = funcml_train,
             model = spec$model, seed = 20260427, trace = FALSE),
        if (!is.null(spec$extra)) spec$extra else list()
      )
      do.call(funcml::fit, fit_args)
    }, error = function(e2) e2)
  })

  if (inherits(fitted, "error")) {
    fit_status[[nm]] <- data.frame(learner = nm, display = DISPLAY[nm],
                                   model = spec$model, fitted = FALSE,
                                   note = conditionMessage(fitted), stringsAsFactors = FALSE)
    next
  }

  pp <- tryCatch({
    cp  <- predict_prob(fitted, funcml_cal)
    rp  <- predict_prob(fitted, funcml_test)
    cal_fit <- fit_calibrator(calibration$label, cp)
    list(cal = cp, raw = rp, pred = apply_calibrator(cal_fit, rp), calibrator = cal_fit)
  }, error = function(e) e)

  if (inherits(pp, "error")) {
    fit_status[[nm]] <- data.frame(learner = nm, display = DISPLAY[nm],
                                   model = spec$model, fitted = FALSE,
                                   note = conditionMessage(pp), stringsAsFactors = FALSE)
    next
  }

  models[[nm]]       <- fitted; funcml_models[[nm]] <- fitted
  calib_preds[[nm]]  <- pp$cal; raw_preds[[nm]] <- pp$raw
  preds[[nm]]        <- pp$pred; calibrators[[nm]] <- pp$calibrator
  fit_status[[nm]]   <- data.frame(learner = nm, display = DISPLAY[nm],
                                   model = spec$model, fitted = TRUE,
                                   note = "ok", stringsAsFactors = FALSE)
  message("  Done: AUROC=", round(auroc(y_test, pp$pred), 3),
          " AUPRC=", round(auprc(y_test, pp$pred), 3))
}

if (length(preds) < 4) stop("Fewer than four learners fitted successfully.", call. = FALSE)

# ── Performance table ─────────────────────────────────────────────────────────
performance <- do.call(rbind, lapply(names(preds), function(nm)
  compute_metrics(y_test, preds[[nm]], nm)))
performance <- performance[order(-performance$AUPRC), ]

# ── Best hyperparameters table ────────────────────────────────────────────────
best_params_rows <- lapply(names(best_params_list), function(nm) {
  bp <- best_params_list[[nm]]
  if (is.null(bp) || nrow(bp) == 0) return(NULL)
  cbind(learner = nm, display = DISPLAY[nm], bp)
})
best_params_df <- do.call(rbind, Filter(Negate(is.null), best_params_rows))

# ── Imputation audit ──────────────────────────────────────────────────────────
imp_audit <- data.frame(
  predictor        = predictors,
  n_missing_train  = vapply(raw_train[, predictors, drop=FALSE], function(x) sum(is.na(x)), integer(1)),
  n_missing_cal    = vapply(raw_cal[,   predictors, drop=FALSE], function(x) sum(is.na(x)), integer(1)),
  n_missing_test   = vapply(raw_test[,  predictors, drop=FALSE], function(x) sum(is.na(x)), integer(1)),
  method           = "mimar RF fit on training window only",
  stringsAsFactors = FALSE
)

# ── Prediction table for downstream scripts ───────────────────────────────────
pred_table <- data.frame(variant_id      = test$variant_id,
                         label           = y_test,
                         final_class_group = test$final_class_group)
for (nm in names(preds))     pred_table[[nm]]                  <- safe_prob(preds[[nm]])
for (nm in names(raw_preds)) pred_table[[paste0(nm, "_raw")]]  <- safe_prob(raw_preds[[nm]])

# ── Save everything ───────────────────────────────────────────────────────────
dir.create("outputs/models",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/predictions", recursive = TRUE, showWarnings = FALSE)

model_bundle <- list(
  models        = models,
  funcml_models = funcml_models,
  tuning_results = tuning_results,
  calibrators   = calibrators,
  predictors    = predictors,
  funcml_levels = funcml_levels,
  display_names = DISPLAY
)
saveRDS(model_bundle,  "outputs/models/all_models.rds")
saveRDS(tuning_results,"outputs/models/tuning_results.rds")
saveRDS(list(models        = imp_models,
             numeric_preds = numeric_preds,
             cols_to_imp   = cols_to_imp,
             apply_fn      = apply_imputations),
        "outputs/models/mimar_rf_imputer.rds")

utils::write.csv(performance,   "outputs/tables/model_performance.csv",      row.names = FALSE)
utils::write.csv(do.call(rbind, fit_status),
                                "outputs/tables/learner_fit_status.csv",      row.names = FALSE)
if (!is.null(best_params_df) && nrow(best_params_df) > 0)
  utils::write.csv(best_params_df, "outputs/tables/best_hyperparameters.csv", row.names = FALSE)
utils::write.csv(imp_audit,     "outputs/tables/imputation_audit.csv",       row.names = FALSE)
utils::write.csv(pred_table,    "outputs/predictions/temporal_test_predictions_all_learners.csv",
                                                                              row.names = FALSE)

message("[06] Done: ", length(preds), " learners fitted. Best AUPRC=",
        round(max(performance$AUPRC, na.rm=TRUE), 3),
        " (", performance$display[1], ").")
