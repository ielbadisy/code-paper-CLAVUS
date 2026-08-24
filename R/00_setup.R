# Package versions used to produce the reported results:
#   R            4.4.1
#   funcml       0.9.0   (CRAN)
#   mimar        1.0.0   (CRAN)
#   glmnet       5.0
#   ranger       0.18.0
#   xgboost      3.2.1.1
#   gbm          2.3.1
#   earth        5.3.5
#   e1071        1.7.17
#   nnet         7.3.20
#   kknn         1.4.1
#   dbarts       0.9.33
#   data.table   1.18.4
#   ggplot2      4.0.3
#   patchwork    1.3.2
#   ragg         1.4.0
#   knitr        1.51
#   quarto       1.5.1

set.seed(20260427)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

required_dirs <- c(
  "data/raw/clinvar", "data/raw/dbnsfp", "data/raw/external", "data/interim",
  "data/curated", "outputs/tables", "outputs/figures", "outputs/models",
  "outputs/predictions", "outputs/logs", "reports", "tests"
)
invisible(lapply(file.path(project_root, required_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

log_message <- function(path, message, append = TRUE) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message),
      file = file.path(project_root, path), append = append)
}

load_or_stop <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this pipeline.", package), call. = FALSE)
  }
}

optional_pkg <- function(package) requireNamespace(package, quietly = TRUE)

clinical_classes <- list(
  p_lp = c("Pathogenic", "Likely pathogenic"),
  b_lb = c("Benign", "Likely benign"),
  excluded = c("Uncertain significance", "Conflicting interpretations of pathogenicity")
)

class_group <- function(x) {
  ifelse(x %in% clinical_classes$p_lp, "P_LP",
         ifelse(x %in% clinical_classes$b_lb, "B_LB", NA_character_))
}

review_stars <- function(x) {
  out <- rep(0, length(x))
  out[grepl("criteria provided, single submitter", x, ignore.case = TRUE)] <- 1
  out[grepl("criteria provided, multiple submitters, no conflicts", x, ignore.case = TRUE)] <- 2
  out[grepl("expert panel", x, ignore.case = TRUE)] <- 3
  out[grepl("practice guideline", x, ignore.case = TRUE)] <- 4
  out
}

safe_prob <- function(x) pmin(pmax(as.numeric(x), 1e-6), 1 - 1e-6)

write.csv2 <- utils::write.csv

log_message("outputs/logs/session_info.txt", capture.output(sessionInfo()), append = FALSE)
log_message("outputs/logs/pipeline_warnings.txt", "Pipeline started. Real ClinVar archived variant_summary releases are required; no synthetic data fallback is used.", append = FALSE)
