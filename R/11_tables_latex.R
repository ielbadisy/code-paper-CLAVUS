# ==============================================================================
# 11_tables_latex.R  —  Write LaTeX table fragments from results CSVs
#
# Produces .tex files in outputs/tables/tex/ for \input{} in the manuscript.
# Each table uses: booktabs + threeparttable (notes) + siunitx (number alignment)
# No kable(); all LaTeX is hand-written for full typographic control.
# ==============================================================================

source("R/00_setup.R")

TEX_DIR <- "outputs/tables/tex"
TBL_DIR <- "outputs/tables"
dir.create(TEX_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ───────────────────────────────────────────────────────────────────
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  sub("(\\w)", "\\U\\1", x, perl = TRUE)
}

LEARNER_DISPLAY <- c(
  xgboost = "XGBoost", bart = "BART", ranger = "Random Forest", gbm = "GBM",
  earth = "MARS", glmnet = "Elastic Net", svm = "Radial SVM",
  nnet = "Neural Network", knn = "k-NN", super_learner = "Super Learner"
)
disp_learner <- function(x) ifelse(x %in% names(LEARNER_DISPLAY), LEARNER_DISPLAY[x], clean_name(x))

fmt3 <- function(x) sprintf("%.3f", as.numeric(x))
fmtN <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)

# Write a complete .tex table file
write_tex <- function(filename, body, caption, label, notes = NULL) {
  header <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\begin{threeparttable}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}")
  )
  footer <- c(
    if (!is.null(notes)) c("\\begin{tablenotes}\\small", paste0("\\item ", notes), "\\end{tablenotes}"),
    "\\end{threeparttable}",
    "\\end{table}"
  )
  writeLines(c(header, body, footer),
             file.path(TEX_DIR, filename))
  message("[tex] Wrote: ", filename)
}

# ── Table 1 — Cohort Characteristics ─────────────────────────────────────────
dat <- readRDS("data/curated/vus_resolution_dataset.rds")
n_total <- nrow(dat); n_plp <- sum(dat$label==1); n_blb <- sum(dat$label==0)

cohort_rows <- c(
  "\\begin{tabular}{lS[table-format=6.0]S[table-format=4.0]S[table-format=4.0]}",
  "\\toprule",
  "Characteristic & {Total} & {P/LP resolved} & {B/LB resolved} \\\\",
  "\\midrule",
  paste0("Variants & ", fmtN(n_total), " & ", fmtN(n_plp), " & ", fmtN(n_blb), " \\\\"),
  paste0("Unique genes & ", fmtN(length(unique(dat$gene_symbol))),
         " & {---} & {---} \\\\"),
  paste0("Median days to resolution & ",
         round(median(dat$days_since_first_seen, na.rm=TRUE)),
         " & ", round(median(dat$days_since_first_seen[dat$label==1], na.rm=TRUE)),
         " & ", round(median(dat$days_since_first_seen[dat$label==0], na.rm=TRUE)), " \\\\"),
  paste0("Star rating {$\\geq$}2 (\\%) & ",
         round(100*mean(dat$baseline_stars_numeric>=2, na.rm=TRUE),1),
         " & ", round(100*mean(dat$baseline_stars_numeric[dat$label==1]>=2, na.rm=TRUE),1),
         " & ", round(100*mean(dat$baseline_stars_numeric[dat$label==0]>=2, na.rm=TRUE),1),
         " \\\\"),
  paste0("gnomAD AF available (\\%) & ",
         round(100*mean(!is.na(dat$gnomad_af)),1), " & {---} & {---} \\\\"),
  paste0("CADD score available (\\%) & ",
         round(100*mean(!is.na(dat$cadd_score)),1), " & {---} & {---} \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex("table1_cohort_characteristics.tex", cohort_rows,
          caption = "Cohort characteristics for the curated longitudinal VUS reclassification dataset.",
          label   = "tab:cohort",
          notes   = paste0("n = ", fmtN(n_total), " variants from six archived ClinVar releases (Jan 2020 -- Apr 2026). ",
                           "P/LP = Pathogenic or Likely Pathogenic; B/LB = Benign or Likely Benign. ",
                           "Days to resolution measured from first ClinVar submission to resolution date."))

# ── Table 2 — Learner Comparison ─────────────────────────────────────────────
perf <- read.csv(file.path(TBL_DIR, "model_performance.csv"),
                 stringsAsFactors = FALSE)
perf <- perf[perf$learner != "super_learner", ]
perf <- perf[order(-perf$AUPRC), ]
perf$display <- disp_learner(perf$learner)

lc_rows <- c(
  "\\begin{tabular}{lS[table-format=1.3]S[table-format=1.3]S[table-format=1.3]S[table-format=1.3]}",
  "\\toprule",
  "Learner & {AUROC} & {AUPRC} & {Sensitivity} & {Specificity} \\\\",
  "\\midrule",
  sapply(seq_len(nrow(perf)), function(i) {
    r <- perf[i, ]
    bold <- r$learner == readLines("outputs/models/best_model_name.txt")[1]
    nm   <- if (bold) paste0("\\textbf{", r$display, "}") else r$display
    paste0(nm, " & ", fmt3(r$AUROC), " & ", fmt3(r$AUPRC), " & ",
           fmt3(r$sensitivity), " & ", fmt3(r$specificity), " \\\\")
  }),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex("table2_learner_comparison.tex", lc_rows,
          caption = "Performance of nine supervised learning algorithms on the temporal test set, ordered by AUPRC.",
          label   = "tab:learner_comparison",
          notes   = "Temporal test set (chronological 70/30 split). Best learner in bold. Threshold = 0.50 for Sensitivity and Specificity. All learners tuned by time-aware rolling cross-validation; metrics computed on the held-out partition without refitting.")

# ── Table 3 — Validation Schemes ─────────────────────────────────────────────
val  <- tryCatch(read.csv(file.path(TBL_DIR, "best_model_validation.csv"),
                           stringsAsFactors=FALSE),
                 error = function(e) read.csv(file.path(TBL_DIR, "table4_best_model_validation.csv"),
                                               stringsAsFactors=FALSE))
val$Scheme <- c("Temporal test", "Gene held-out")[seq_len(nrow(val))]

best_nm_tbl <- readLines("outputs/models/best_model_name.txt")[1]
best_disp   <- if (best_nm_tbl %in% names(LEARNER_DISPLAY)) LEARNER_DISPLAY[[best_nm_tbl]] else clean_name(best_nm_tbl)
boot_m_tbl  <- read.csv(file.path(TBL_DIR, "bootstrapped_metrics_best_model.csv"), stringsAsFactors=FALSE)
auroc_lo    <- sprintf("%.3f", boot_m_tbl$ci_lower[boot_m_tbl$metric == "AUROC"])
auroc_hi    <- sprintf("%.3f", boot_m_tbl$ci_upper[boot_m_tbl$metric == "AUROC"])
auroc_ci_str <- paste0(auroc_lo, "--", auroc_hi)

vl_rows <- c(
  "\\begin{tabular}{lS[table-format=5.0]lS[table-format=1.3]S[table-format=1.3]S[table-format=1.3]S[table-format=1.3]}",
  "\\toprule",
  "Validation scheme & {n} & {AUROC (95\\% CI)} & {AUPRC} & {Brier} & {Sensitivity} & {Specificity} \\\\",
  "\\midrule",
  sapply(seq_len(nrow(val)), function(i) {
    r <- val[i, ]
    auroc_cell <- if (i == 1)
      paste0(fmt3(r$AUROC), " (", auroc_ci_str, ")")
    else
      fmt3(r$AUROC)
    paste0(r$Scheme, " & ", fmtN(r$n), " & ", auroc_cell, " & ",
           fmt3(r$AUPRC), " & ", fmt3(r$Brier_score), " & ",
           fmt3(r$sensitivity), " & ", fmt3(r$specificity), " \\\\")
  }),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex("table3_validation.tex", vl_rows,
          caption = paste0(best_disp, " model performance under temporal and gene-held-out validation schemes."),
          label   = "tab:validation",
          notes   = paste0("Temporal: chronological 70/30 split; Gene held-out: all variants from 20\\% of genes withheld from training and model selection (seed 20260427). AUROC bootstrapped 95\\% CI (B = 1000): temporal ", auroc_lo, "--", auroc_hi, ", gene held-out not reported."))

# ── Table 4 — External Annotation Coverage ────────────────────────────────────
cov <- read.csv(file.path(TBL_DIR, "external_annotation_coverage.csv"),
                stringsAsFactors=FALSE)
cov$annotation <- clean_name(cov$annotation)

cv_rows <- c(
  "\\begin{tabular}{lS[table-format=5.0]S[table-format=2.1]}",
  "\\toprule",
  "Annotation & {n available} & {\\% available} \\\\",
  "\\midrule",
  sapply(seq_len(nrow(cov)), function(i) {
    r <- cov[i,]
    paste0(r$annotation, " & ", fmtN(r$n_available), " & ",
           sprintf("%.1f", r$pct_available), " \\\\")
  }),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex("table4_annotation_coverage.tex", cv_rows,
          caption = "External functional annotation coverage for the curated VUS cohort.",
          label   = "tab:coverage",
          notes   = paste0("Total cohort n = ", fmtN(nrow(dat)),
                           ". Coverage = non-missing values returned by MyVariant.info for GRCh38 coordinates."))

# ── Table 5 — ATE Results ─────────────────────────────────────────────────────
ate <- read.csv(file.path(TBL_DIR, "ate_results.csv"), stringsAsFactors=FALSE)
METHOD_DISPLAY <- c(ipw="IPW", aipw="AIPW (doubly robust)", tmle="TMLE (one-step)")
exposures_ord <- unique(ate$exposure)

ate_rows <- c(
  paste0("\\begin{tabular}{p{4.5cm}",
         paste(rep("S[table-format=+1.3]", 3), collapse=""),
         "l}"),
  "\\toprule",
  paste0("Exposure & {IPW} & {AIPW} & {TMLE} & {Agreement} \\\\"),
  "\\midrule"
)
for (exp in exposures_ord) {
  sub <- ate[ate$exposure == exp, ]
  agr <- sub$sign_agreement[1]
  vals <- sapply(c("ipw","aipw","tmle"), function(m) {
    r <- sub[sub$method==m, ]
    if (nrow(r)==0 || is.na(r$estimate)) return("---")
    paste0(sprintf("%.3f",r$estimate))
  })
  agr_cell <- if (agr=="agree") "\\checkmark" else "\\ding{55}"
  ate_rows <- c(ate_rows,
    paste0("\\multicolumn{5}{l}{\\textit{", exp, "}} \\\\"),
    paste0(" & ", paste(vals, collapse=" & "), " & ", agr_cell, " \\\\")
  )
}
ate_rows <- c(ate_rows, "\\bottomrule", "\\end{tabular}")

n_disagree  <- sum(unique(ate[, c("exposure","sign_agreement")])$sign_agreement == "disagree")
write_tex("table5_ate_results.tex", ate_rows,
          caption = "Adjusted-effect estimates (risk differences) for four binary exposures on P/LP resolution probability. Estimates shown for three estimators (IPW, AIPW, TMLE); CIs from Rubin-pooled bootstrap (Track B, m = 5) or bootstrap (Track A, B = 200).",
          label   = "tab:ate",
          notes   = paste0("IPW = inverse probability weighting (stabilized, truncated); AIPW = augmented IPW (doubly robust); TMLE = targeted maximum likelihood. \\checkmark = all three estimators agree on direction of effect; \\ding{55} = sign disagreement across estimators (", n_disagree, " of 4 exposures), indicating sensitivity to unmeasured confounding."))

# ── Table 6 — Best Hyperparameters ───────────────────────────────────────────
if (file.exists(file.path(TBL_DIR, "best_hyperparameters.csv"))) {
  hp <- read.csv(file.path(TBL_DIR, "best_hyperparameters.csv"), stringsAsFactors=FALSE)
  hp_param_cols <- setdiff(names(hp), c("learner","display"))
  hp_rows <- c(
    paste0("\\begin{tabular}{l", paste(rep("l", length(hp_param_cols)), collapse=""), "}"),
    "\\toprule",
    paste0("Learner & ", paste(clean_name(hp_param_cols), collapse=" & "), " \\\\"),
    "\\midrule",
    sapply(seq_len(nrow(hp)), function(i) {
      r <- hp[i,]
      vals <- sapply(hp_param_cols, function(col) {
        v <- r[[col]]
        if (is.na(v)) "---" else as.character(v)
      })
      paste0(r$display, " & ", paste(vals, collapse=" & "), " \\\\")
    }),
    "\\bottomrule",
    "\\end{tabular}"
  )
  write_tex("table6_best_hyperparameters.tex", hp_rows,
            caption = "Optimal hyperparameters selected by 3-fold cross-validated Brier score minimisation for each learner.",
            label   = "tab:hyperparameters",
            notes   = "All learners tuned on the training partition (n = 70\\% of cohort). Tuning grid sizes: XGBoost 16, BART 3, Random Forest 8, GBM 8, MARS 4, Elastic Net 9, Radial SVM 6, Neural Network 4, k-NN 3.")
}

message("[11_tables] All LaTeX table fragments written to: ", TEX_DIR)
