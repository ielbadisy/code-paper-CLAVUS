# ==============================================================================
# 11_figures_production.R  —  Nature Medicine Figures (400 DPI, ragg + cairo)
#
# All 8 figures. Standards:
#   • ragg::agg_png  400 DPI, width 183 mm (Nature double-column)
#   • cairo_pdf      vector backup
#   • No underscores in any axis label, legend text, or title
#   • Okabe-Ito colorblind-safe palette, consistent mapping
#   • theme_classic(base_size = 10, base_family = "sans")
#   • patchwork multi-panel, panel labels A/B/C
# ==============================================================================

source("R/00_setup.R")
library(ggplot2)
library(patchwork)
library(ragg)

FIG_DIR  <- "outputs/figures"
TBL_DIR  <- "outputs/tables"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Okabe-Ito palette — named ─────────────────────────────────────────────────
OI <- c(
  orange    = "#E69F00", sky_blue  = "#56B4E9", green     = "#009E73",
  yellow    = "#F0E442", blue      = "#0072B2", vermillion = "#D55E00",
  pink      = "#CC79A7", black     = "#000000", grey      = "#999999"
)

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_nm <- function(...) {
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(face = "bold", size = 10),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      ...
    )
}

# ── Save helpers ──────────────────────────────────────────────────────────────
save_fig <- function(p, name, width = 183, height = 120, units = "mm") {
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  ragg::agg_png(png_path, width = width, height = height, units = units, res = 400)
  print(p); dev.off()
  cairo_pdf(pdf_path, width = width / 25.4, height = height / 25.4)
  print(p); dev.off()
  message("[fig] Saved: ", name)
}

# ── Display-name helpers ──────────────────────────────────────────────────────
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)  # Title Case
  x
}

LEARNER_DISPLAY <- c(
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
disp <- function(x) {
  ifelse(x %in% names(LEARNER_DISPLAY), LEARNER_DISPLAY[x], clean_name(x))
}

# ==============================================================================
# Figure 1 — Study Design & Cohort Flow
# ==============================================================================
dat       <- readRDS("data/curated/vus_resolution_dataset.rds")
n_total   <- nrow(dat)
n_plp     <- sum(dat$final_class_group == "P_LP")
n_blb     <- sum(dat$final_class_group == "B_LB")
n_genes   <- length(unique(dat$gene_symbol))

releases  <- c("Jan 2020\n(baseline)", "Jan 2021", "Jan 2022",
                "Jan 2023", "Dec 2024", "Apr 2026")
x_pos     <- seq_len(length(releases))
y_rel     <- rep(2, length(releases))

rel_df <- data.frame(x = x_pos, y = y_rel, label = releases)
arr_df <- data.frame(x    = x_pos[-length(x_pos)] + 0.05,
                     xend = x_pos[-1]              - 0.05,
                     y = 2, yend = 2)

cohort_text <- data.frame(
  x = c(1.4, 3.5, 3.5),
  y = c(0.7, 1.3, 0.3),
  label = c(
    paste0("Baseline VUS cohort\n(n = ", format(n_total, big.mark=","), " variants\nacross ", n_genes, " genes)"),
    paste0("Resolved P/LP\nn = ", format(n_plp, big.mark=","), " (", round(100*n_plp/n_total,1), "%)"),
    paste0("Resolved B/LB\nn = ", format(n_blb, big.mark=","), " (", round(100*n_blb/n_total,1), "%)")
  ),
  col = c(OI["grey"], OI["vermillion"], OI["sky_blue"])
)

p_flow <- ggplot() +
  geom_segment(data = arr_df, aes(x=x, xend=xend, y=y, yend=yend),
               arrow = arrow(length = unit(0.2,"cm"), type="closed"),
               linewidth = 0.7, colour = OI["black"]) +
  geom_point(data = rel_df, aes(x=x, y=y), size=5,
             colour = OI["blue"], shape = 21, fill = OI["sky_blue"], stroke=1.2) +
  geom_text(data = rel_df, aes(x=x, y=y+0.18, label=label),
            size=2.8, fontface="bold", hjust=0.5) +
  geom_label(data = cohort_text, aes(x=x, y=y, label=label),
             fill = "white", colour = cohort_text$col,
             size = 2.5, label.size = 0.4, label.r = unit(0.1,"cm")) +
  scale_y_continuous(limits = c(0, 2.5)) +
  labs(title = "A  ClinVar archival releases and VUS reclassification cohort",
       x = "ClinVar Release", y = NULL) +
  theme_nm() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), panel.grid = element_blank())

# Prevalence bar inset
prev_df <- data.frame(
  Class     = c("Pathogenic / Likely Pathogenic", "Benign / Likely Benign"),
  Proportion = c(n_plp, n_blb) / n_total,
  Fill = c(OI["vermillion"], OI["sky_blue"])
)
p_bar <- ggplot(prev_df, aes(x="", y=Proportion, fill=Class)) +
  geom_col(width=0.5, colour="white") +
  geom_text(aes(label=paste0(round(Proportion*100,1),"%")),
            position=position_stack(vjust=0.5), size=3, colour="white", fontface="bold") +
  scale_fill_manual(values=setNames(prev_df$Fill, prev_df$Class)) +
  labs(title="B  Resolution class\ndistribution", x=NULL, y="Proportion", fill=NULL) +
  coord_flip() + theme_nm() +
  theme(legend.position="bottom", axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

fig1 <- p_flow / p_bar + plot_layout(heights=c(3,1))
save_fig(fig1, "figure1_study_design_cohort", height=140)

# ==============================================================================
# Figure 2 — Learner Comparison (AUROC + AUPRC dot plot)
# ==============================================================================
perf <- read.csv(file.path(TBL_DIR, "model_performance.csv"),
                 stringsAsFactors = FALSE)
perf <- perf[perf$learner != "super_learner", ]
perf$display <- disp(perf$learner)
perf$display <- factor(perf$display, levels = rev(perf$display[order(perf$AUPRC)]))

melt_perf <- rbind(
  data.frame(display=perf$display, Metric="AUROC", Value=perf$AUROC),
  data.frame(display=perf$display, Metric="AUPRC", Value=perf$AUPRC)
)

p_cmp <- ggplot(melt_perf, aes(x=Value, y=display, colour=Metric, shape=Metric)) +
  geom_vline(xintercept=0.5, linetype="dashed", colour="grey70", linewidth=0.4) +
  geom_point(size=3.5) +
  scale_colour_manual(values=c(AUROC=OI[["blue"]], AUPRC=OI[["vermillion"]])) +
  scale_shape_manual(values=c(AUROC=16, AUPRC=17)) +
  scale_x_continuous(limits=c(0.45, 1), breaks=seq(0.5,1,0.1)) +
  labs(x="Metric value", y=NULL, colour=NULL, shape=NULL,
       title="Learner comparison on temporal test set") +
  theme_nm() +
  theme(legend.position=c(0.85, 0.15))

save_fig(p_cmp, "figure2_learner_comparison", height=100)

# ==============================================================================
# Figure 3 — ROC and PR Curves (best model)
# ==============================================================================
pred_tbl <- read.csv("outputs/predictions/temporal_test_predictions_all_learners.csv",
                     stringsAsFactors = FALSE)
best_nm  <- readLines("outputs/models/best_model_name.txt")[1]
y        <- pred_tbl$label
p_best   <- pred_tbl[[best_nm]]

# Bootstrap CIs (authoritative values)
boot_m <- read.csv(file.path(TBL_DIR, "bootstrapped_metrics_best_model.csv"),
                   stringsAsFactors = FALSE)
fmt3 <- function(x) sprintf("%.3f", x)
# Point estimates from direct test-set evaluation; CIs from bootstrap
val_tbl     <- read.csv(file.path(TBL_DIR, "best_model_validation.csv"), stringsAsFactors=FALSE)
val_tbl     <- val_tbl[val_tbl$validation_scheme == "temporal_test", ]
auroc_val   <- fmt3(val_tbl$AUROC)
auroc_ci_lo <- fmt3(boot_m$ci_lower[boot_m$metric == "AUROC"])
auroc_ci_hi <- fmt3(boot_m$ci_upper[boot_m$metric == "AUROC"])
auprc_val   <- fmt3(val_tbl$AUPRC)
auprc_ci_lo <- fmt3(boot_m$ci_lower[boot_m$metric == "AUPRC"])
auprc_ci_hi <- fmt3(boot_m$ci_upper[boot_m$metric == "AUPRC"])

# ROC
thresholds <- sort(unique(c(0, p_best, 1)))
roc_df <- do.call(rbind, lapply(thresholds, function(t) {
  pred_bin <- as.integer(p_best >= t)
  tp <- sum(pred_bin == 1 & y == 1); fp <- sum(pred_bin == 1 & y == 0)
  tn <- sum(pred_bin == 0 & y == 0); fn <- sum(pred_bin == 0 & y == 1)
  data.frame(threshold=t, TPR=tp/max(tp+fn,1), FPR=fp/max(fp+tn,1))
}))
roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]

p_roc <- ggplot(roc_df, aes(x=FPR, y=TPR)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60", linewidth=0.5) +
  geom_line(colour=OI["blue"], linewidth=0.9) +
  annotate("text", x=0.55, y=0.15,
           label=paste0("AUROC = ", auroc_val,
                        "\n(95% CI ", auroc_ci_lo, "–", auroc_ci_hi, ")"),
           size=3.2, hjust=0, colour=OI["blue"]) +
  scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)",
       title="A  Receiver Operating Characteristic") +
  theme_nm()

# PR curve
ord <- order(p_best, decreasing=TRUE)
y_ord  <- y[ord]; p_ord <- p_best[ord]
tp_cum <- cumsum(y_ord == 1); fp_cum <- cumsum(y_ord == 0)
pr_df  <- data.frame(
  Recall    = tp_cum / max(sum(y==1),1),
  Precision = tp_cum / pmax(tp_cum + fp_cum, 1)
)
null_prec <- mean(y)

p_pr <- ggplot(pr_df, aes(x=Recall, y=Precision)) +
  geom_hline(yintercept=null_prec, linetype="dashed", colour="grey60", linewidth=0.5) +
  geom_line(colour=OI["vermillion"], linewidth=0.9) +
  annotate("text", x=0.05, y=0.15,
           label=paste0("AUPRC = ", auprc_val,
                        "\n(95% CI ", auprc_ci_lo, "–", auprc_ci_hi, ")"),
           size=3.2, hjust=0, colour=OI["vermillion"]) +
  scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  labs(x="Recall (Sensitivity)", y="Precision (PPV)",
       title="B  Precision–Recall Curve") +
  theme_nm()

fig3 <- p_roc + p_pr
save_fig(fig3, "figure3_roc_pr_curves", height=90)

# ==============================================================================
# Figure 4 — Calibration Reliability Diagram (best model)
# ==============================================================================
calib_diag <- read.csv(file.path(TBL_DIR,"calibration_diagnostics.csv"),
                        stringsAsFactors=FALSE)
# Use selected_best_model flag — learner IDs in this CSV may differ from best_nm
gb_diag    <- calib_diag[isTRUE(calib_diag$selected_best_model) |
                          calib_diag$selected_best_model == "TRUE", ][1, ]

n_bins <- 10
cuts   <- cut(p_best, breaks=seq(0,1,length.out=n_bins+1), include.lowest=TRUE)
cal_df <- aggregate(cbind(obs=y, pred=p_best), list(bin=cuts), mean)
cal_df$n <- as.numeric(table(cuts)[as.character(cal_df$bin)])

loess_fit <- stats::loess(obs ~ pred, data=cal_df, weights=cal_df$n, span=1)
pred_seq  <- seq(0.01, 0.99, length.out=100)
loess_df  <- data.frame(pred=pred_seq, obs=predict(loess_fit, pred_seq))
loess_df  <- loess_df[!is.na(loess_df$obs),]

ece_v <- if(nrow(gb_diag)>0) sprintf("%.3f", gb_diag$ECE) else "—"
slp_v <- if(nrow(gb_diag)>0) sprintf("%.3f", gb_diag$calibration_slope) else "—"

p_cal <- ggplot() +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=0.5) +
  geom_line(data=loess_df, aes(x=pred, y=obs), colour=OI["blue"], linewidth=1) +
  geom_point(data=cal_df, aes(x=pred, y=obs, size=n), colour=OI["vermillion"],
             alpha=0.85, shape=21, fill=OI["orange"], stroke=0.8) +
  scale_size_area(max_size=6, guide="none") +
  annotate("text", x=0.65, y=0.10,
           label=paste0("ECE = ", ece_v, "\nCalibration slope = ", slp_v),
           size=3.2, hjust=0, colour=OI["blue"]) +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  labs(x="Mean predicted probability", y="Observed event rate",
       title=paste0("Calibration — ", disp(best_nm), " (after Platt recalibration)")) +
  theme_nm()

save_fig(p_cal, "figure4_calibration_curve_best_model", width=150, height=110)

# ==============================================================================
# Figure 5 — Variable Importance
# ==============================================================================
vip <- read.csv(file.path(TBL_DIR,"table5_top_predictors.csv"), stringsAsFactors=FALSE)
vip$display_var <- clean_name(vip$variable)
vip$display_var <- factor(vip$display_var, levels=rev(vip$display_var[order(vip$importance)]))

top_n <- min(20, nrow(vip))
vip   <- vip[seq_len(top_n), ]

p_vip <- ggplot(vip, aes(x=importance, y=display_var)) +
  geom_col(fill=OI["blue"], alpha=0.85, width=0.7) +
  scale_x_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x="Variable importance (gain)", y=NULL,
       title=paste0("Variable importance — ", disp(best_nm), " (top ", top_n, " features)")) +
  theme_nm() +
  theme(panel.grid.major.y=element_blank())

save_fig(p_vip, "figure5_variable_importance_best_model", height=130)

# ==============================================================================
# Figure 6 — Decision Curve Analysis
# ==============================================================================
if (file.exists(file.path(TBL_DIR, "dca_results.csv"))) {
  dca_raw <- read.csv(file.path(TBL_DIR, "dca_results.csv"), stringsAsFactors=FALSE)

  # Reshape: expect wide (outcome, threshold, <model>, treat_all, treat_none)
  if (!"strategy" %in% names(dca_raw)) {
    model_cols <- setdiff(names(dca_raw), c("outcome","threshold","treat_all","treat_none"))
    dca_long <- do.call(rbind, c(
      lapply(model_cols, function(mc)
        data.frame(outcome=dca_raw$outcome, threshold=dca_raw$threshold,
                   strategy=disp(mc), net_benefit=dca_raw[[mc]],
                   type="Model", stringsAsFactors=FALSE)),
      list(data.frame(outcome=dca_raw$outcome, threshold=dca_raw$threshold,
                      strategy="Treat All",  net_benefit=dca_raw$treat_all, type="Reference",stringsAsFactors=FALSE)),
      list(data.frame(outcome=dca_raw$outcome, threshold=dca_raw$threshold,
                      strategy="Treat None", net_benefit=dca_raw$treat_none,type="Reference",stringsAsFactors=FALSE))
    ))
  } else {
    dca_long <- dca_raw
    dca_long$strategy <- clean_name(dca_long$strategy)
  }

  outcomes_present <- unique(dca_long$outcome)
  dca_list <- lapply(outcomes_present, function(oc) {
    sub <- dca_long[dca_long$outcome == oc & !is.na(dca_long$net_benefit), ]
    sub <- sub[sub$threshold >= 0 & sub$threshold <= 0.95, ]
    refs <- c("Treat All","Treat None")
    sub$linetype_val <- ifelse(sub$strategy %in% refs, sub$strategy, "Model")

    ggplot(sub, aes(x=threshold, y=net_benefit,
                    colour=strategy, linetype=linetype_val, linewidth=linetype_val)) +
      geom_line() +
      scale_colour_manual(values=c(
        setNames(OI["blue"], unique(sub$strategy[!sub$strategy %in% refs])[1]),
        "Treat All"  = OI[["orange"]],
        "Treat None" = OI[["grey"]]
      ), na.value=OI["grey"]) +
      scale_linetype_manual(values=c(Model="solid","Treat All"="longdash","Treat None"="dotted"),
                             guide="none") +
      scale_linewidth_manual(values=c(Model=1.0,"Treat All"=0.7,"Treat None"=0.7),
                              guide="none") +
      coord_cartesian(ylim = c(
        -0.3,
        max(sub$net_benefit[sub$strategy != "Treat All"], na.rm=TRUE) + 0.02
      )) +
      scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
      labs(x="Decision threshold", y="Net benefit", colour=NULL,
           title=paste0(if(grepl("P_LP|plp",oc,ignore.case=TRUE)) "A  P/LP prioritization"
                        else "B  B/LB triage")) +
      theme_nm() +
      theme(legend.position=c(0.75,0.75))
  })
  fig6 <- Reduce(`+`, dca_list) + plot_layout(ncol=2)
  save_fig(fig6, "figure6_decision_curve_analysis", height=100)
}

# ==============================================================================
# Figure 7 — Predicted Probability Distribution by Final Class
# ==============================================================================
score_df <- data.frame(
  Score = p_best,
  Class = factor(ifelse(y==1, "P/LP resolved", "B/LB resolved"),
                 levels = c("P/LP resolved", "B/LB resolved"))
)
p_dist <- ggplot(score_df, aes(x=Score, fill=Class, colour=Class)) +
  geom_density(alpha=0.35, linewidth=0.8) +
  geom_vline(xintercept=c(0.20, 0.80), linetype="dashed",
             colour=OI["black"], linewidth=0.6) +
  scale_fill_manual(values=c("P/LP resolved"=OI[["vermillion"]], "B/LB resolved"=OI[["sky_blue"]])) +
  scale_colour_manual(values=c("P/LP resolved"=OI[["vermillion"]], "B/LB resolved"=OI[["sky_blue"]])) +
  annotate("text", x=c(0.10,0.50,0.90), y=Inf, vjust=1.5, size=2.8,
           label=c("Likely B/LB","Uncertain","Likely P/LP"), colour=OI["black"]) +
  labs(x="Predicted P/LP resolution probability", y="Density",
       fill=NULL, colour=NULL,
       title="Predicted score distribution by resolved classification") +
  theme_nm() +
  theme(legend.position=c(0.85,0.85))

save_fig(p_dist, "figure7_prediction_score_distribution_by_final_class", width=120, height=100)

# ==============================================================================
# Figure 8 — ATE Forest Plot (rebuilt, high quality)
# ==============================================================================
ate <- read.csv(file.path(TBL_DIR,"ate_results.csv"), stringsAsFactors=FALSE)

METHOD_DISPLAY <- c(
  ipw  = "IPW",
  aipw = "AIPW",
  tmle = "TMLE"
)
METHOD_COLORS <- c(
  "IPW"  = OI[["orange"]],
  "AIPW" = OI[["green"]],
  "TMLE" = OI[["pink"]]
)
ate$method_display <- METHOD_DISPLAY[ate$method]

# Clip extreme CIs (unstable IPW bootstrap) to ±1 risk difference
ate$conf_low  <- pmax(ate$conf_low,  ate$estimate - 1)
ate$conf_high <- pmin(ate$conf_high, ate$estimate + 1)

# Agreed exposures (panel A) vs disagreed (panel B)
agreed    <- ate[ate$sign_agreement == "agree",    ]
disagreed <- ate[ate$sign_agreement == "disagree", ]

make_forest_panel <- function(df, title, alpha_val=1, show_legend=TRUE) {
  if (nrow(df) == 0) return(NULL)
  # levels must match the display VALUES (not the key names)
  df$method_display <- factor(df$method_display, levels=rev(unname(METHOD_DISPLAY)))
  df$exposure_short <- gsub("evidence|probability|above median at t0|at the|non-","",df$exposure)
  df$exposure_short <- trimws(gsub("  +"," ",df$exposure_short))
  df$exposure_short <- factor(df$exposure_short, levels=unique(df$exposure_short))

  ggplot(df, aes(x=estimate, y=method_display, colour=method_display,
                  xmin=conf_low, xmax=conf_high)) +
    geom_vline(xintercept=0, colour="grey40", linewidth=0.5) +
    geom_errorbarh(height=0.2, alpha=alpha_val, linewidth=0.8) +
    geom_point(size=3, alpha=alpha_val) +
    facet_wrap(~exposure_short, ncol=1, scales="free_x") +
    scale_colour_manual(values=METHOD_COLORS, guide=if(show_legend) "legend" else "none") +
    labs(x="Average Treatment Effect (risk difference)", y=NULL, colour="Estimator",
         title=title) +
    theme_nm() +
    theme(legend.position=if(show_legend) "right" else "none",
          strip.text=element_text(size=8.5, face="bold"),
          panel.spacing=unit(0.6,"lines"))
}

p_agreed    <- make_forest_panel(agreed,    "A  Agreed-direction exposures (IPW, AIPW, TMLE)")
p_disagreed <- make_forest_panel(disagreed, "B  Sign-disagreement exposures\n(not interpreted causally)",
                                  alpha_val=0.45, show_legend=FALSE)

if (!is.null(p_agreed) && !is.null(p_disagreed)) {
  fig8 <- p_agreed + p_disagreed + plot_layout(ncol=2, widths=c(3,2))
} else {
  fig8 <- p_agreed %||% p_disagreed
}
save_fig(fig8, "figure8_ate_forest_plot", height=160)

message("[11_figures] All figures saved.")
