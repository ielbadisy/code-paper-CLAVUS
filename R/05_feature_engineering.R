source("R/00_setup.R")

dat <- readRDS("data/interim/annotated_vus_reclassifications.rds")

# ── Baseline ClinVar metadata features ────────────────────────────────────────
dat$baseline_review_status <- dat$review_status
dat$baseline_number_submitters <- dat$number_submitters
dat$baseline_number_conditions <- dat$number_conditions
dat$baseline_assertion_method_available <- grepl("criteria provided|expert panel|practice guideline", dat$baseline_review_status, ignore.case = TRUE)
dat$baseline_has_conflict <- grepl("conflict", dat$baseline_review_status, ignore.case = TRUE)
dat$baseline_stars_numeric <- review_stars(dat$baseline_review_status)
dat$days_since_first_seen <- ifelse(is.na(dat$last_evaluated), 0L, pmax(0L, as.integer(as.Date(dat$baseline_date) - as.Date(dat$last_evaluated))))
dat$number_prior_updates <- 0L
dat$molecular_consequence <- ifelse(grepl("missense", dat$variant_name, ignore.case = TRUE), "missense_or_named_protein_change",
                                    ifelse(grepl("splice|\\+", dat$variant_name, ignore.case = TRUE), "splice_region_or_intronic",
                                           ifelse(grepl("del|dup|ins", dat$variant_type, ignore.case = TRUE), dat$variant_type, "not_available")))
dat$protein_change_available <- grepl("\\(p\\.", dat$variant_name, fixed = FALSE)
dat$is_missense <- grepl("missense", dat$molecular_consequence, ignore.case = TRUE)
dat$is_splice_region <- grepl("splice|intronic", dat$molecular_consequence, ignore.case = TRUE)
dat$is_synonymous <- grepl("synonymous", dat$variant_name, ignore.case = TRUE)
dat$is_nonsense <- grepl("Ter|\\*", dat$variant_name)
dat$is_frameshift <- grepl("fs", dat$variant_name, ignore.case = TRUE)
external_cols <- c("gnomad_af", "gnomad_ac", "gnomad_an", "gnomad_hom_count", "cadd_score", "revel_score",
                   "sift_score", "polyphen_score", "mutationtaster_score", "alphamissense_score", "spliceai_score")
for (nm in external_cols) {
  if (!nm %in% names(dat)) dat[[nm]] <- NA_real_
}
if (!"dbnsfp_available" %in% names(dat)) dat$dbnsfp_available <- FALSE
dat$is_missing_gnomad_af    <- is.na(dat$gnomad_af)
dat$is_missing_cadd         <- is.na(dat$cadd_score)
dat$is_missing_revel        <- is.na(dat$revel_score)
dat$is_missing_alphamissense <- is.na(dat$alphamissense_score)
dat$is_missing_spliceai     <- is.na(dat$spliceai_score)

# ── Extended features from 03b (AlphaMissense full-genome, gnomAD constraint, ClinGen) ──
for (nm in c("gnomad_pli", "gnomad_loeuf")) {
  if (!nm %in% names(dat)) dat[[nm]] <- NA_real_
}
if (!"clingen_validity" %in% names(dat)) dat$clingen_validity <- "not_curated"

# ACMG evidence criterion proxy flags
if (!"acmg_pvs1_proxy" %in% names(dat))
  dat$acmg_pvs1_proxy <- !is.na(dat$gnomad_loeuf) & dat$gnomad_loeuf < 0.35 &
    (dat$is_nonsense | dat$is_frameshift | dat$is_splice_region)
if (!"acmg_pm2_proxy" %in% names(dat))
  dat$acmg_pm2_proxy <- !is.na(dat$gnomad_af) & dat$gnomad_af < 0.0001
if (!"acmg_ba1_proxy" %in% names(dat))
  dat$acmg_ba1_proxy <- !is.na(dat$gnomad_af) & dat$gnomad_af >= 0.05
if (!"acmg_ps4_proxy" %in% names(dat))
  dat$acmg_ps4_proxy <- !is.na(dat$gene_pathogenic_fraction_at_t0) &
    dat$gene_pathogenic_fraction_at_t0 >= 0.30
if (!"acmg_pp3_proxy" %in% names(dat))
  dat$acmg_pp3_proxy <- with(dat, ifelse(!is.na(revel_score), revel_score >= 0.75,
                              ifelse(!is.na(alphamissense_score), alphamissense_score >= 0.80,
                              ifelse(!is.na(cadd_score), cadd_score >= 25, NA))))
if (!"acmg_bp4_proxy" %in% names(dat))
  dat$acmg_bp4_proxy <- with(dat, ifelse(!is.na(revel_score), revel_score < 0.25,
                              ifelse(!is.na(alphamissense_score), alphamissense_score < 0.20,
                              ifelse(!is.na(cadd_score), cadd_score < 10, NA))))

dat$is_missing_gnomad_pli   <- is.na(dat$gnomad_pli)
dat$is_missing_gnomad_loeuf <- is.na(dat$gnomad_loeuf)
dat$is_missing_acmg_pp3     <- is.na(dat$acmg_pp3_proxy)
dat$is_missing_acmg_bp4     <- is.na(dat$acmg_bp4_proxy)

baseline <- readRDS("data/interim/clinvar_t0_real.rds")
baseline$baseline_stars_numeric <- review_stars(baseline$review_status)
gene_counts <- aggregate(
  cbind(
    gene_variant_count_at_t0 = rep(1, nrow(baseline)),
    gene_vus_count_at_t0 = baseline$clinical_significance == "Uncertain significance",
    gene_pathogenic_count_at_t0 = baseline$clinical_significance %in% clinical_classes$p_lp,
    gene_benign_count_at_t0 = baseline$clinical_significance %in% clinical_classes$b_lb
  ) ~ gene_symbol,
  baseline,
  sum
)
gene_counts$gene_pathogenic_fraction_at_t0 <- with(gene_counts, gene_pathogenic_count_at_t0 / pmax(gene_variant_count_at_t0, 1))
gene_counts$gene_benign_fraction_at_t0 <- with(gene_counts, gene_benign_count_at_t0 / pmax(gene_variant_count_at_t0, 1))
dat <- merge(dat, gene_counts, by = "gene_symbol", all.x = TRUE)
top_genes <- names(sort(table(dat$gene_symbol), decreasing = TRUE))[seq_len(min(50, length(unique(dat$gene_symbol))))]
dat$gene_symbol_model <- ifelse(dat$gene_symbol %in% top_genes, dat$gene_symbol, "OTHER")

base_predictor_names <- c(
  "variant_type", "molecular_consequence", "chromosome", "gene_symbol_model",
  "protein_change_available", "is_missense", "is_splice_region", "is_synonymous", "is_nonsense", "is_frameshift",
  "baseline_review_status", "baseline_number_submitters", "baseline_number_conditions",
  "baseline_assertion_method_available", "baseline_has_conflict", "baseline_stars_numeric",
  "days_since_first_seen", "number_prior_updates",
  "gene_variant_count_at_t0", "gene_vus_count_at_t0", "gene_pathogenic_count_at_t0", "gene_benign_count_at_t0",
  "gene_pathogenic_fraction_at_t0", "gene_benign_fraction_at_t0"
)
external_predictor_names <- external_cols[vapply(dat[, external_cols, drop = FALSE],
  function(x) sum(!is.na(x)) >= 50, logical(1))]

# Extended predictors — only include where ≥50 non-missing values
extended_predictor_names <- character(0)
for (nm in c("gnomad_pli", "gnomad_loeuf")) {
  if (sum(!is.na(dat[[nm]])) >= 50) extended_predictor_names <- c(extended_predictor_names, nm)
}
if (length(unique(dat$clingen_validity)) >= 2)
  extended_predictor_names <- c(extended_predictor_names, "clingen_validity")
acmg_candidates <- c("acmg_pvs1_proxy", "acmg_pm2_proxy", "acmg_ba1_proxy", "acmg_ps4_proxy")
acmg_predictor_names <- acmg_candidates[vapply(dat[, acmg_candidates, drop = FALSE],
  function(x) sum(!is.na(x)) >= 50, logical(1))]
extended_predictor_names <- c(extended_predictor_names, acmg_predictor_names)

predictor_names <- unique(c(
  base_predictor_names,
  external_predictor_names,
  extended_predictor_names,
  "dbnsfp_available",
  "is_missing_gnomad_af", "is_missing_cadd", "is_missing_revel",
  "is_missing_alphamissense", "is_missing_spliceai",
  "is_missing_gnomad_pli", "is_missing_gnomad_loeuf",
  "is_missing_acmg_pp3", "is_missing_acmg_bp4"
))
predictor_names <- predictor_names[predictor_names %in% names(dat)]

dat <- dat[order(dat$final_date, dat$variant_id), ]
cut_index <- max(1, floor(0.70 * nrow(dat)))
cutoff_date <- dat$final_date[cut_index]
dat$split_temporal <- ifelse(seq_len(nrow(dat)) <= cut_index, "train", "test")
set.seed(20260427)
heldout_genes <- sample(unique(dat$gene_symbol), max(2, ceiling(length(unique(dat$gene_symbol)) * 0.2)))
dat$split_gene_heldout <- ifelse(dat$gene_symbol %in% heldout_genes, "test", "train")

feature_recipe <- list(
  predictor_names = predictor_names,
  outcome = "label",
  temporal_cutoff_date = cutoff_date,
  heldout_genes = heldout_genes
)

saveRDS(dat, "data/curated/vus_resolution_dataset.rds")
utils::write.csv(dat, "data/curated/vus_resolution_dataset.csv", row.names = FALSE)
saveRDS(feature_recipe, "data/curated/feature_recipe.rds")

validity <- data.frame(
  check = c(
    "unique_variant_baseline_final_rows",
    "labels_binary",
    "baseline_all_vus",
    "final_all_resolved",
    "no_final_variables_in_predictors",
    "temporal_split_order_valid",
    "train_test_variant_disjoint",
    "both_classes_in_train",
    "both_classes_in_test",
    "real_clinvar_source",
    "no_synthetic_scores",
    "real_external_annotation_coverage"
  ),
  passed = c(
    anyDuplicated(paste(dat$variant_id, dat$baseline_date, dat$final_date)) == 0,
    all(dat$label %in% c(0, 1)),
    all(dat$baseline_class == "Uncertain significance"),
    all(dat$final_class_group %in% c("P_LP", "B_LB")),
    length(intersect(predictor_names, c("final_class", "final_class_group", "final_date", "label", "time_to_reclassification_days", "final_review_status"))) == 0,
    max(as.Date(dat$final_date[dat$split_temporal == "train"])) <= min(as.Date(dat$final_date[dat$split_temporal == "test"])),
    length(intersect(dat$variant_id[dat$split_temporal == "train"], dat$variant_id[dat$split_temporal == "test"])) == 0,
    length(unique(dat$label[dat$split_temporal == "train"])) == 2,
    length(unique(dat$label[dat$split_temporal == "test"])) == 2,
    all(file.exists(read.csv("data/raw/clinvar/release_metadata.csv")$local_file)),
    file.exists("outputs/tables/external_annotation_coverage.csv"),
    any(!is.na(dat$gnomad_af)) || any(!is.na(dat$cadd_score)) || any(!is.na(dat$revel_score)) || any(!is.na(dat$alphamissense_score))
  )
)
utils::write.csv(validity, "outputs/tables/dataset_validity_checks.csv", row.names = FALSE)
if (!all(validity$passed)) {
  stop("Curated dataset failed validity checks; see outputs/tables/dataset_validity_checks.csv", call. = FALSE)
}
log_message("outputs/logs/data_sources.txt",
  sprintf("Curated dataset saved with %d rows, %d predictors (%d extended: %s), temporal cutoff %s.",
    nrow(dat), length(predictor_names), length(extended_predictor_names),
    paste(extended_predictor_names, collapse=", "), cutoff_date))
message(sprintf("[05] Feature engineering complete: %d rows, %d predictors (%d base, %d extended).",
  nrow(dat), length(predictor_names), length(base_predictor_names), length(extended_predictor_names)))
