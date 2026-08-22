# ==============================================================================
# 10_ate_estimation.R  —  ATE Estimation with Two-Track Imputation Strategy
#
# Track A (full-cohort, ClinVar-native covariates — no external annotations):
#   Complete-case; no imputation needed.
#
# Track B (exposures relying on external annotations — CADD, gnomAD):
#   m=20 mimar RF imputations via mimar::impute(); pooling via mimar::pool().
#
# Four estimators per exposure:
#   1. G-computation      (funcml::estimate outcome model)
#   2. IPW                (stabilized, 1st/99th-pctile truncated, bootstrap SE)
#   3. AIPW               (doubly robust, bootstrap SE)
#   4. TMLE               (one-step, efficient-IF SE)
#
# Outputs:
#   outputs/tables/ate_results.csv
#   outputs/tables/ate_sign_agreement.csv
#   outputs/tables/ate_evalue_summary.csv
# ==============================================================================

source("R/00_setup.R")
library(funcml)
library(mimar)

# ── Load models and data ──────────────────────────────────────────────────────
dat        <- readRDS("data/curated/vus_resolution_dataset.rds")
model_bndl <- readRDS("outputs/models/all_models.rds")

best_model_name <- readLines("outputs/models/best_model_name.txt")[1]
best_model      <- model_bndl$funcml_models[[best_model_name]]
best_calibrator <- model_bndl$calibrators[[best_model_name]]
predictors      <- model_bndl$predictors
funcml_levels   <- model_bndl$funcml_levels

BASELINE_RISK <- mean(dat$label)
M_IMPUTATIONS <- 5L
N_BOOT        <- 200L

logit  <- function(p) log(p / (1 - p))
expit  <- function(x) 1 / (1 + exp(-x))
safe_p <- function(p) pmax(pmin(as.numeric(p), 1 - 1e-6), 1e-6)

# ── Exposure definitions (applied to any data frame) ─────────────────────────
define_exposures <- function(d) {
  d$high_deleteriousness <- with(d, ifelse(
    !is.na(revel_score),            revel_score >= 0.75,
    ifelse(!is.na(alphamissense_score), alphamissense_score >= 0.80,
      ifelse(!is.na(cadd_score),        cadd_score >= 25, NA))))
  d$nonrare_frequency   <- ifelse(is.na(d$gnomad_af), NA, d$gnomad_af >= 0.001)
  d$high_review_support <- d$baseline_stars_numeric >= 2
  d$gene_pathogenic_enriched <- d$gene_pathogenic_fraction_at_t0 >=
    stats::median(d$gene_pathogenic_fraction_at_t0, na.rm = TRUE)
  d
}
dat <- define_exposures(dat)

# ── Exposure registry ─────────────────────────────────────────────────────────
EXPOSURES <- list(
  list(var = "high_deleteriousness",     label = "High computational deleteriousness evidence",
       expected_sign = "positive",  track = "B"),
  list(var = "nonrare_frequency",        label = "Population-frequency evidence non-rare",
       expected_sign = "negative",  track = "B"),
  list(var = "high_review_support",      label = "High baseline ClinVar review support",
       expected_sign = NA_character_, track = "A"),
  list(var = "gene_pathogenic_enriched", label = "Gene pathogenic-enrichment above median at t0",
       expected_sign = NA_character_, track = "A")
)

BASE_COVARS <- c(
  "baseline_number_submitters", "baseline_number_conditions",
  "gene_vus_count_at_t0",       "gene_benign_fraction_at_t0",
  "variant_type",               "molecular_consequence",
  "chromosome",
  "is_missing_gnomad_af",       "is_missing_cadd",
  "is_missing_revel",           "is_missing_alphamissense"
)

# External numeric columns that mimar will impute for Track B
EXT_COLS <- c("gnomad_af", "gnomad_ac", "gnomad_an", "gnomad_hom_count",
              "cadd_score", "revel_score", "alphamissense_score",
              "sift_score", "polyphen_score")
EXT_COLS <- EXT_COLS[EXT_COLS %in% names(dat)]

# ── mimar Track-B: generate m=20 imputed datasets ────────────────────────────
# Only columns relevant to ATE are passed to mimar to keep memory manageable
ATE_COLS <- unique(c("label", "high_deleteriousness", "nonrare_frequency",
                     BASE_COVARS, EXT_COLS))
ATE_COLS <- ATE_COLS[ATE_COLS %in% names(dat)]

message("[10] Generating m=", M_IMPUTATIONS, " mimar RF imputations ...")
mi_obj     <- mimar::impute(dat[, ATE_COLS], m = M_IMPUTATIONS,
                             imputer = "rf", seed = 20260427, verbose = FALSE)
mi_datasets <- lapply(seq_len(M_IMPUTATIONS), function(m) {
  d_m <- as.data.frame(mimar::complete(mi_obj, m))
  # Re-attach columns not in ATE_COLS needed for prediction
  extra <- setdiff(names(dat), ATE_COLS)
  cbind(d_m, dat[, extra, drop = FALSE])
})
# Re-define exposures on each imputed dataset
mi_datasets <- lapply(mi_datasets, define_exposures)
message("[10] MI datasets ready.")

# ── Outcome prediction helper ─────────────────────────────────────────────────
prep_funcml <- function(d) {
  out <- d[, predictors, drop = FALSE]
  for (nm in predictors) {
    if (is.logical(out[[nm]])) out[[nm]] <- as.integer(out[[nm]])
    if (is.character(out[[nm]])) out[[nm]][is.na(out[[nm]])] <- "Missing"
    if (nm %in% names(funcml_levels)) {
      lv <- funcml_levels[[nm]]
      v  <- as.character(out[[nm]]); v[!v %in% lv] <- "OTHER"
      out[[nm]] <- factor(v, levels = lv)
    }
    if (is.numeric(out[[nm]]) && any(is.na(out[[nm]]))) {
      med <- stats::median(out[[nm]], na.rm = TRUE)
      out[[nm]][is.na(out[[nm]])] <- if (is.na(med)) 0 else med
    }
  }
  out
}

pred_outcome <- function(d) {
  safe_p(stats::predict(best_calibrator,
    newdata = data.frame(p = safe_p(
      stats::predict(best_model, newdata = prep_funcml(d), type = "prob")[, "P_LP"])),
    type = "response"))
}

# ── Core estimators ───────────────────────────────────────────────────────────
fit_ps <- function(d, evar, covars) {
  A <- as.integer(d[[evar]])
  ps <- if (length(covars) > 0) {
    f  <- stats::as.formula(paste(evar, "~", paste(covars, collapse = "+")))
    m  <- tryCatch(stats::glm(f, data = d, family = stats::binomial()),
                   error = function(e) NULL)
    if (!is.null(m)) safe_p(stats::predict(m, type = "response")) else rep(mean(A), nrow(d))
  } else rep(mean(A), nrow(d))
  p_m <- mean(A)
  w   <- ifelse(A == 1L, p_m / ps, (1 - p_m) / (1 - ps))
  q   <- stats::quantile(w, c(0.01, 0.99))
  list(ps = ps, w = pmin(pmax(w, q[1]), q[2]), A = A)
}

est_gcomp <- function(d, evar) {
  d1 <- d; d1[[evar]] <- 1L
  d0 <- d; d0[[evar]] <- 0L
  mean(pred_outcome(d1) - pred_outcome(d0))
}

est_ipw <- function(d, evar, covars, y, B = N_BOOT) {
  ps_obj <- fit_ps(d, evar, covars); A <- ps_obj$A; w <- ps_obj$w
  ate    <- mean(A * y * w) - mean((1 - A) * y * w)
  boots  <- vapply(seq_len(B), function(b) {
    set.seed(b); idx <- sample(nrow(d), replace = TRUE)
    db <- d[idx,]; yb <- y[idx]; pso <- fit_ps(db, evar, covars)
    mean(pso$A * yb * pso$w) - mean((1 - pso$A) * yb * pso$w)
  }, numeric(1))
  list(estimate = ate, SE = stats::sd(boots))
}

est_aipw <- function(d, evar, covars, y, B = N_BOOT) {
  ps_obj <- fit_ps(d, evar, covars); A <- ps_obj$A; e_x <- ps_obj$ps
  mu  <- pred_outcome(d)
  d1  <- d; d1[[evar]] <- 1L; d0 <- d; d0[[evar]] <- 0L
  m1  <- pred_outcome(d1); m0 <- pred_outcome(d0)
  dr  <- (m1 - m0) + A * (y - mu) / e_x - (1 - A) * (y - mu) / (1 - e_x)
  ate <- mean(dr)
  boots <- vapply(seq_len(B), function(b) {
    set.seed(b + 1000L); idx <- sample(nrow(d), replace = TRUE)
    db <- d[idx,]; yb <- y[idx]; pso <- fit_ps(db, evar, covars)
    Ab <- pso$A; eb <- pso$ps; mub <- pred_outcome(db)
    db1 <- db; db1[[evar]] <- 1L; db0 <- db; db0[[evar]] <- 0L
    m1b <- pred_outcome(db1); m0b <- pred_outcome(db0)
    mean((m1b - m0b) + Ab * (yb - mub) / eb - (1 - Ab) * (yb - mub) / (1 - eb))
  }, numeric(1))
  list(estimate = ate, SE = stats::sd(boots))
}

est_tmle <- function(d, evar, covars, y) {
  tryCatch({
    ps_obj <- fit_ps(d, evar, covars); A <- ps_obj$A; e_x <- ps_obj$ps
    Q  <- pred_outcome(d); H1 <- A / e_x; H0 <- -(1 - A) / (1 - e_x)
    ep <- tryCatch(unname(coef(stats::glm(y ~ -1 + I(H1 + H0),
                                           offset = logit(Q),
                                           family  = stats::binomial()))[1]),
                   error = function(e) 0)
    d1 <- d; d1[[evar]] <- 1L; d0 <- d; d0[[evar]] <- 0L
    Q1 <- expit(logit(pred_outcome(d1)) + ep / e_x)
    Q0 <- expit(logit(pred_outcome(d0)) - ep / (1 - e_x))
    ate <- mean(Q1 - Q0)
    ic  <- (Q1 - Q0) - ate + A * (y - expit(logit(Q) + ep * H1)) / e_x -
           (1 - A) * (y - expit(logit(Q) + ep * H0)) / (1 - e_x)
    list(estimate = ate, SE = sqrt(stats::var(ic) / nrow(d)))
  }, error = function(e) {
    message("  TMLE failed: ", e$message); list(estimate = NA_real_, SE = NA_real_)
  })
}

# ── E-value ───────────────────────────────────────────────────────────────────
compute_evalue <- function(rd) {
  if (length(rd) == 0 || is.na(rd)) return(NA_real_)
  rr <- (BASELINE_RISK + rd) / BASELINE_RISK
  rr <- max(rr, 1 / rr, na.rm = TRUE)
  rr + sqrt(rr * (rr - 1))
}

# ── Run one exposure on one dataset; return named estimates + SEs ─────────────
run_one <- function(d, evar, B = N_BOOT) {
  cols_ok <- BASE_COVARS[BASE_COVARS %in% names(d)]
  cols_ok <- setdiff(cols_ok, evar)
  sub     <- d[!is.na(d[[evar]]), , drop = FALSE]
  sub     <- sub[stats::complete.cases(sub[, c("label", evar, cols_ok)]), , drop = FALSE]
  covars  <- cols_ok[vapply(sub[, cols_ok, drop=FALSE],
                            function(x) length(unique(x)) >= 2L, logical(1))]
  y       <- sub$label
  list(
    n             = nrow(sub),
    g_computation = list(estimate = tryCatch(est_gcomp(sub, evar), error=function(e) NA_real_),
                         SE = NA_real_),
    ipw           = tryCatch(est_ipw( sub, evar, covars, y, B), error=function(e) list(estimate=NA_real_, SE=NA_real_)),
    aipw          = tryCatch(est_aipw(sub, evar, covars, y, B), error=function(e) list(estimate=NA_real_, SE=NA_real_)),
    tmle          = tryCatch(est_tmle(sub, evar, covars, y),    error=function(e) list(estimate=NA_real_, SE=NA_real_))
  )
}

# ── Pool m results using mimar::pool() ───────────────────────────────────────
pool_method_results <- function(mi_results, method) {
  ests <- vapply(mi_results, function(r) r[[method]]$estimate, numeric(1))
  ses  <- vapply(mi_results, function(r) r[[method]]$SE,       numeric(1))
  if (all(is.na(ses))) {
    # No variance info — use robust pooling
    p <- mimar::pool(ests, name = method)
  } else {
    p <- mimar::pool(ests, std.error = ses, name = method)
  }
  pdata <- p$pooled
  list(
    estimate  = pdata$estimate[1],
    SE        = pdata$std.error[1],
    conf_low  = pdata$conf.low[1],
    conf_high = pdata$conf.high[1]
  )
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_rows <- list()

for (exp_cfg in EXPOSURES) {
  message("[10] Exposure: ", exp_cfg$label, " (Track ", exp_cfg$track, ")")

  if (exp_cfg$track == "B") {
    # Run on each MI dataset with reduced bootstrap (50) for speed, pool with mimar
    B_mi <- 50L
    mi_res <- lapply(seq_along(mi_datasets), function(m) {
      message("  imputation m=", m, "/", M_IMPUTATIONS)
      run_one(mi_datasets[[m]], exp_cfg$var, B = B_mi)
    })
    n_avg  <- round(mean(vapply(mi_res, `[[`, numeric(1), "n")))
    pooled <- lapply(c("g_computation", "ipw", "aipw", "tmle"),
                     function(meth) pool_method_results(mi_res, meth))
    names(pooled) <- c("g_computation", "ipw", "aipw", "tmle")

    # G-computation: bootstrap CI via full N_BOOT on the primary (observed) complete case
    sub0 <- dat[!is.na(dat[[exp_cfg$var]]), , drop = FALSE]
    covok <- setdiff(BASE_COVARS[BASE_COVARS %in% names(sub0)], exp_cfg$var)
    sub0  <- sub0[stats::complete.cases(sub0[, c("label", exp_cfg$var, covok)]), , drop = FALSE]
    gc_boot <- vapply(seq_len(N_BOOT), function(b) {
      set.seed(b + 9000L); idx <- sample(nrow(sub0), replace = TRUE)
      tryCatch(est_gcomp(sub0[idx, ], exp_cfg$var), error = function(e) NA_real_)
    }, numeric(1))
    gc_se  <- stats::sd(gc_boot, na.rm = TRUE)
    gc_est <- pooled$g_computation$estimate
    pooled$g_computation <- list(
      estimate  = gc_est,
      SE        = gc_se,
      conf_low  = gc_est - 1.96 * gc_se,
      conf_high = gc_est + 1.96 * gc_se
    )

  } else {  # Track A — single complete-case run, full bootstrap
    res    <- run_one(dat, exp_cfg$var, B = N_BOOT)
    n_avg  <- res$n
    pooled <- list(
      g_computation = res$g_computation,
      ipw           = res$ipw,
      aipw          = res$aipw,
      tmle          = res$tmle
    )
    # Bootstrap CI for g-computation
    sub0  <- dat[!is.na(dat[[exp_cfg$var]]), , drop = FALSE]
    covok <- setdiff(BASE_COVARS[BASE_COVARS %in% names(sub0)], exp_cfg$var)
    sub0  <- sub0[stats::complete.cases(sub0[, c("label", exp_cfg$var, covok)]), , drop = FALSE]
    gc_boot <- vapply(seq_len(N_BOOT), function(b) {
      set.seed(b + 9000L); idx <- sample(nrow(sub0), replace = TRUE)
      tryCatch(est_gcomp(sub0[idx, ], exp_cfg$var), error = function(e) NA_real_)
    }, numeric(1))
    gc_se  <- stats::sd(gc_boot, na.rm = TRUE)
    gc_est <- pooled$g_computation$estimate
    pooled$g_computation <- list(
      estimate  = gc_est,
      SE        = gc_se,
      conf_low  = gc_est - 1.96 * gc_se,
      conf_high = gc_est + 1.96 * gc_se
    )
  }

  for (meth in c("ipw", "aipw", "tmle")) {
    r   <- pooled[[meth]]
    est <- if (length(r$estimate) > 0) r$estimate else NA_real_
    se  <- if (length(r$SE) > 0) r$SE else NA_real_
    cl  <- if (!is.null(r$conf_low)  && length(r$conf_low)  > 0) r$conf_low  else est - 1.96 * se
    ch  <- if (!is.null(r$conf_high) && length(r$conf_high) > 0) r$conf_high else est + 1.96 * se
    ev  <- compute_evalue(est)
    all_rows[[length(all_rows) + 1]] <- data.frame(
      exposure        = exp_cfg$label,
      method          = meth,
      track           = exp_cfg$track,
      estimate        = est,
      SE              = se,
      conf_low        = cl,
      conf_high       = ch,
      n_complete      = n_avg,
      expected_sign   = exp_cfg$expected_sign,
      evalue          = ev,
      evalue_lower_ci = compute_evalue(cl),
      stringsAsFactors = FALSE
    )
  }
}

ate <- do.call(rbind, all_rows)

# ── Sign agreement ────────────────────────────────────────────────────────────
sign_agree <- do.call(rbind, lapply(unique(ate$exposure), function(exp) {
  sub  <- ate[ate$exposure == exp & !is.na(ate$estimate), ]
  agr  <- length(unique(sign(sub$estimate))) == 1L
  by_m <- paste(paste0(sub$method, "=", round(sub$estimate, 4)), collapse = "; ")
  data.frame(exposure           = exp,
             n_estimators_run   = nrow(sub),
             sign_agreement     = ifelse(agr, "agree", "disagree"),
             estimates_by_method = by_m,
             stringsAsFactors   = FALSE)
}))
ate <- merge(ate, sign_agree[, c("exposure", "sign_agreement")], by = "exposure", all.x = TRUE)

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(ate,        "outputs/tables/ate_results.csv",        row.names = FALSE)
utils::write.csv(sign_agree, "outputs/tables/ate_sign_agreement.csv", row.names = FALSE)
utils::write.csv(ate[ate$method == "aipw",
                     c("exposure", "estimate", "evalue", "evalue_lower_ci", "sign_agreement")],
                 "outputs/tables/ate_evalue_summary.csv", row.names = FALSE)

message("[10] Done. Sign agreement summary:")
print(sign_agree[, c("exposure", "sign_agreement", "estimates_by_method")])
