source("R/00_setup.R")

# Downloads and integrates four additional feature layers:
#   1. AlphaMissense (Google DeepMind, Science 2023) - pre-computed genome-wide TSV
#   2. gnomAD v2.1.1 gene constraint table (LOEUF, pLI)
#   3. ClinGen gene-disease validity curations
#   4. ACMG evidence criterion flags from existing ClinVar submission metadata

message("[03b] Downloading and integrating external feature layers...")

dat <- readRDS("data/interim/annotated_vus_reclassifications.rds")

# ── 1. AlphaMissense ─────────────────────────────────────────────────────────
# Pre-computed genome-wide scores (Cheng et al., Science 2023)
# URL: https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz
# ~1.4 GB — skipped automatically if absent; fallback uses alphamissense_score
# from MyVariant.info already in the dataset

am_path <- "data/raw/external/AlphaMissense_hg38.tsv.gz"
if (!file.exists(am_path)) {
  am_url <- "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"
  message("[03b] AlphaMissense file absent. Attempting download (~1.4 GB)...")
  tryCatch({
    utils::download.file(am_url, am_path, mode = "wb", quiet = FALSE, method = "libcurl")
    message("[03b] AlphaMissense downloaded successfully.")
  }, error = function(e) {
    message("[03b] AlphaMissense download failed: ", conditionMessage(e),
            ". Continuing with MyVariant.info alphamissense_score where available.")
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] AlphaMissense download failed:", conditionMessage(e)))
  })
}

if (file.exists(am_path)) {
  message("[03b] Integrating AlphaMissense scores...")
  tryCatch({
    # Build CHROM:POS:REF:ALT key from dataset
    dat$am_key <- with(dat, paste0(chromosome, ":", position, ":", reference_allele, ":", alternative_allele))
    # Stream AlphaMissense TSV in chunks to avoid memory overflow
    con <- gzcon(file(am_path, "rb"))
    on.exit(close(con), add = TRUE)
    header <- readLines(con, n = 3)
    # Find the actual header line (starts with #CHROM)
    hdr_line <- header[grepl("^#CHROM|^CHROM|^chrom", header, ignore.case = TRUE)]
    if (length(hdr_line) == 0) {
      message("[03b] AlphaMissense header not found in first 3 lines; skipping integration.")
    } else {
      col_names <- strsplit(gsub("^#", "", hdr_line), "\t")[[1]]
      am_scores <- list()
      chunk_size <- 500000L
      repeat {
        lines <- readLines(con, n = chunk_size)
        if (length(lines) == 0) break
        chunk <- data.frame(
          do.call(rbind, strsplit(lines, "\t")),
          stringsAsFactors = FALSE
        )
        names(chunk) <- col_names[seq_len(ncol(chunk))]
        # Standardise column names
        names(chunk) <- tolower(gsub("^#", "", names(chunk)))
        key_col  <- intersect(c("chrom", "chromosome"), names(chunk))[1]
        pos_col  <- intersect(c("pos", "position"), names(chunk))[1]
        ref_col  <- intersect(c("ref", "reference_allele"), names(chunk))[1]
        alt_col  <- intersect(c("alt", "alternative_allele"), names(chunk))[1]
        score_col <- intersect(c("am_pathogenicity", "pathogenicity", "score"), names(chunk))[1]
        if (any(is.na(c(key_col, pos_col, ref_col, alt_col, score_col)))) break
        chunk$am_key_chunk <- paste0(chunk[[key_col]], ":", chunk[[pos_col]], ":",
                                      chunk[[ref_col]], ":", chunk[[alt_col]])
        hits <- chunk[chunk$am_key_chunk %in% dat$am_key,
                      c("am_key_chunk", score_col), drop = FALSE]
        if (nrow(hits) > 0) {
          names(hits) <- c("am_key", "alphamissense_full_score")
          am_scores[[length(am_scores) + 1]] <- hits
        }
      }
      if (length(am_scores) > 0) {
        am_all <- unique(do.call(rbind, am_scores))
        am_all$alphamissense_full_score <- suppressWarnings(as.numeric(am_all$alphamissense_full_score))
        dat <- merge(dat, am_all, by = "am_key", all.x = TRUE)
        # Prefer full-genome score over MyVariant.info partial score
        dat$alphamissense_score <- ifelse(
          !is.na(dat$alphamissense_full_score),
          dat$alphamissense_full_score,
          dat$alphamissense_score
        )
        n_am <- sum(!is.na(dat$alphamissense_full_score))
        message(sprintf("[03b] AlphaMissense: %d/%d variants annotated (%.1f%%)",
                        n_am, nrow(dat), 100 * n_am / nrow(dat)))
        log_message("outputs/logs/annotation_warnings.txt",
                    sprintf("[03b] AlphaMissense full-genome integration: %d/%d variants annotated.", n_am, nrow(dat)))
      }
      dat$am_key <- NULL
      if ("alphamissense_full_score" %in% names(dat)) dat$alphamissense_full_score <- NULL
    }
  }, error = function(e) {
    message("[03b] AlphaMissense integration error: ", conditionMessage(e))
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] AlphaMissense integration error:", conditionMessage(e)))
  })
} else {
  message("[03b] AlphaMissense file not available; using MyVariant.info alphamissense_score where present.")
}

# ── 2. gnomAD v2.1.1 gene constraint (LOEUF, pLI) ───────────────────────────
# Public download: https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
gnomad_constraint_path <- "data/raw/external/gnomad_gene_constraint.tsv"
if (!file.exists(gnomad_constraint_path)) {
  gnomad_url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
  message("[03b] gnomAD gene constraint file absent. Attempting download...")
  tryCatch({
    tmp_gz <- paste0(gnomad_constraint_path, ".bgz")
    utils::download.file(gnomad_url, tmp_gz, mode = "wb", quiet = FALSE, method = "libcurl")
    # bgz is gzip-compatible
    gz_con <- gzcon(file(tmp_gz, "rb"))
    lines <- readLines(gz_con); close(gz_con)
    writeLines(lines, gnomad_constraint_path)
    unlink(tmp_gz)
    message("[03b] gnomAD gene constraint downloaded and extracted.")
  }, error = function(e) {
    message("[03b] gnomAD constraint download failed: ", conditionMessage(e), ". Skipping LOEUF/pLI.")
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] gnomAD gene constraint download failed:", conditionMessage(e)))
  })
}

if (file.exists(gnomad_constraint_path)) {
  message("[03b] Integrating gnomAD gene constraint (LOEUF, pLI)...")
  tryCatch({
    gc_raw <- utils::read.table(gnomad_constraint_path, header = TRUE, sep = "\t",
                                 quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
    gc_cols <- intersect(c("gene", "pLI", "oe_lof_upper", "oe_mis_upper", "oe_syn_upper",
                            "exp_lof", "obs_lof", "oe_lof"), names(gc_raw))
    gc <- gc_raw[, gc_cols, drop = FALSE]
    gc <- gc[!is.na(gc$gene), ]
    names(gc)[names(gc) == "gene"] <- "gene_symbol"
    names(gc)[names(gc) == "pLI"] <- "gnomad_pli"
    names(gc)[names(gc) == "oe_lof_upper"] <- "gnomad_loeuf"
    gc <- gc[!duplicated(gc$gene_symbol), ]
    keep_gc <- c("gene_symbol",
                 intersect(c("gnomad_pli", "gnomad_loeuf", "oe_lof", "oe_mis_upper"), names(gc)))
    dat <- merge(dat, gc[, keep_gc, drop = FALSE], by = "gene_symbol", all.x = TRUE)
    n_loeuf <- sum(!is.na(dat$gnomad_loeuf))
    n_pli   <- sum(!is.na(dat$gnomad_pli))
    message(sprintf("[03b] gnomAD constraint: LOEUF for %d/%d variants, pLI for %d/%d variants.",
                    n_loeuf, nrow(dat), n_pli, nrow(dat)))
  }, error = function(e) {
    message("[03b] gnomAD constraint integration error: ", conditionMessage(e))
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] gnomAD gene constraint integration error:", conditionMessage(e)))
  })
}

# Ensure columns exist even if download failed
for (col in c("gnomad_pli", "gnomad_loeuf")) {
  if (!col %in% names(dat)) dat[[col]] <- NA_real_
}

# ── 3. ClinGen gene-disease validity curations ────────────────────────────────
# Public TSV: https://search.clinicalgenome.org/kb/gene-validity/download
clingen_path <- "data/raw/external/clingen_gene_validity.csv"
if (!file.exists(clingen_path)) {
  clingen_url <- "https://search.clinicalgenome.org/kb/gene-validity/download"
  message("[03b] ClinGen gene validity file absent. Attempting download...")
  tryCatch({
    utils::download.file(clingen_url, clingen_path, mode = "w", quiet = FALSE, method = "libcurl")
    message("[03b] ClinGen gene validity downloaded.")
  }, error = function(e) {
    message("[03b] ClinGen download failed: ", conditionMessage(e), ". Skipping ClinGen validity.")
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] ClinGen download failed:", conditionMessage(e)))
  })
}

if (file.exists(clingen_path)) {
  message("[03b] Integrating ClinGen gene-disease validity...")
  tryCatch({
    cg_raw <- tryCatch(
      utils::read.csv(clingen_path, stringsAsFactors = FALSE, skip = 3),
      error = function(e) utils::read.csv(clingen_path, stringsAsFactors = FALSE)
    )
    # Standardise columns — ClinGen CSV has varying header names across versions
    names(cg_raw) <- tolower(gsub("[ .]", "_", names(cg_raw)))
    gene_col <- intersect(c("gene_symbol", "gene", "hgnc_gene_symbol"), names(cg_raw))[1]
    validity_col <- intersect(c("classification", "validity_classification",
                                 "disease_validity_classification"), names(cg_raw))[1]
    if (!is.na(gene_col) && !is.na(validity_col)) {
      cg <- cg_raw[, c(gene_col, validity_col), drop = FALSE]
      names(cg) <- c("gene_symbol", "clingen_validity_raw")
      # Consolidate to ordered factor levels
      validity_map <- c(
        "Definitive" = "definitive", "Strong" = "strong", "Moderate" = "moderate",
        "Limited" = "limited", "Disputed" = "disputed", "Refuted" = "refuted",
        "Animal Model Only" = "animal_model_only", "No Reported Evidence" = "no_evidence"
      )
      cg$clingen_validity <- validity_map[cg$clingen_validity_raw]
      # Keep strongest classification per gene (Definitive > Strong > Moderate > Limited > ...)
      validity_order <- c("definitive", "strong", "moderate", "limited",
                           "disputed", "animal_model_only", "no_evidence", "refuted")
      cg$validity_rank <- match(cg$clingen_validity, validity_order)
      cg <- cg[order(cg$validity_rank, na.last = TRUE), ]
      cg <- cg[!duplicated(cg$gene_symbol), c("gene_symbol", "clingen_validity")]
      dat <- merge(dat, cg, by = "gene_symbol", all.x = TRUE)
      dat$clingen_validity[is.na(dat$clingen_validity)] <- "not_curated"
      n_cg <- sum(dat$clingen_validity != "not_curated")
      message(sprintf("[03b] ClinGen validity: %d/%d variants in curated genes.", n_cg, nrow(dat)))
    } else {
      message("[03b] ClinGen CSV column names not recognized; skipping.")
    }
  }, error = function(e) {
    message("[03b] ClinGen integration error: ", conditionMessage(e))
    log_message("outputs/logs/annotation_warnings.txt",
                paste("[03b] ClinGen integration error:", conditionMessage(e)))
  })
}

if (!"clingen_validity" %in% names(dat)) dat$clingen_validity <- "not_curated"

# ── 4. ACMG evidence criterion flags from ClinVar metadata ───────────────────
# These are approximated from review_status and submitter count — the full
# ClinVar XML carries per-criterion flags (PVS1, PS1, PM1...) but the
# tab-delimited variant_summary does not. We use the available proxies.

message("[03b] Deriving ACMG evidence proxy flags from ClinVar metadata...")

# PVS1 proxy: loss-of-function consequence in a gene with LOEUF < 0.35
dat$acmg_pvs1_proxy <- !is.na(dat$gnomad_loeuf) &
  dat$gnomad_loeuf < 0.35 &
  (dat$is_nonsense | dat$is_frameshift | dat$is_splice_region)

# PM2 proxy: very low gnomAD allele frequency (rare in population)
dat$acmg_pm2_proxy <- !is.na(dat$gnomad_af) & dat$gnomad_af < 0.0001

# BA1 proxy: allele frequency > 5% — strong benign evidence
dat$acmg_ba1_proxy <- !is.na(dat$gnomad_af) & dat$gnomad_af >= 0.05

# PP3/BP4 proxy: computational evidence — use REVEL/AlphaMissense/CADD
dat$acmg_pp3_proxy <- with(dat, ifelse(!is.na(revel_score), revel_score >= 0.75,
                             ifelse(!is.na(alphamissense_score), alphamissense_score >= 0.80,
                             ifelse(!is.na(cadd_score), cadd_score >= 25, NA))))
dat$acmg_bp4_proxy <- with(dat, ifelse(!is.na(revel_score), revel_score < 0.25,
                             ifelse(!is.na(alphamissense_score), alphamissense_score < 0.20,
                             ifelse(!is.na(cadd_score), cadd_score < 10, NA))))

# PS4 proxy: gene has high pathogenic variant fraction at baseline (established disease gene)
dat$acmg_ps4_proxy <- !is.na(dat$gene_pathogenic_fraction_at_t0) &
  dat$gene_pathogenic_fraction_at_t0 >= 0.30

# Missingness indicators for new features
dat$is_missing_gnomad_pli   <- is.na(dat$gnomad_pli)
dat$is_missing_gnomad_loeuf <- is.na(dat$gnomad_loeuf)
dat$is_missing_acmg_pp3     <- is.na(dat$acmg_pp3_proxy)
dat$is_missing_acmg_bp4     <- is.na(dat$acmg_bp4_proxy)

# ── 5. Save enriched interim dataset ─────────────────────────────────────────
saveRDS(dat, "data/interim/annotated_vus_reclassifications.rds")

# Coverage summary for new features
new_features <- c("alphamissense_score", "gnomad_loeuf", "gnomad_pli",
                   "clingen_validity", "acmg_pvs1_proxy", "acmg_pm2_proxy",
                   "acmg_ba1_proxy", "acmg_pp3_proxy", "acmg_bp4_proxy", "acmg_ps4_proxy")
coverage_new <- data.frame(
  annotation = new_features,
  n_available = vapply(dat[, new_features, drop = FALSE], function(x) {
    if (is.logical(x) || is.numeric(x)) sum(!is.na(x)) else sum(!is.na(x) & x != "not_curated")
  }, integer(1)),
  n_total = nrow(dat),
  stringsAsFactors = FALSE
)
coverage_new$pct_available <- round(100 * coverage_new$n_available / coverage_new$n_total, 2)
utils::write.csv(coverage_new, "outputs/tables/external_feature_coverage_extended.csv", row.names = FALSE)
log_message("outputs/logs/annotation_warnings.txt",
            sprintf("[03b] Extended feature integration complete: %d rows, %d new features.",
                    nrow(dat), length(new_features)))
message("[03b] External feature integration complete.")
