source("R/00_setup.R")
load_or_stop("data.table")

selected_cols <- c(
  "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
  "ClinSigSimple", "LastEvaluated", "RS# (dbSNP)", "Origin", "OriginSimple",
  "Assembly", "Chromosome", "Start", "Stop", "ReferenceAllele",
  "AlternateAllele", "ReviewStatus", "NumberSubmitters", "VariationID"
)

normalize_class <- function(x) {
  x <- trimws(x)
  x <- sub("\\s*/\\s*", "/", x)
  x <- ifelse(x %in% c("Pathogenic/Likely pathogenic", "Pathogenic; Likely pathogenic"), "Pathogenic", x)
  x <- ifelse(x %in% c("Benign/Likely benign", "Benign; Likely benign"), "Benign", x)
  x
}

read_release_dt <- function(path, release_date) {
  x <- data.table::fread(
    path, sep = "\t", quote = "", na.strings = c("", "-", "not provided"),
    select = selected_cols, showProgress = TRUE
  )
  data.table::setnames(x, "#AlleleID", "AlleleID")
  x[, variant_id := sprintf("VCV%09d", as.integer(VariationID))]
  x[, release_date := as.Date(release_date)]
  x[, ClinicalSignificance := normalize_class(ClinicalSignificance)]
  x[, LastEvaluatedDate := suppressWarnings(as.Date(LastEvaluated, format = "%b %d, %Y"))]
  x <- x[
    !is.na(variant_id) &
      !is.na(ClinicalSignificance) &
      (grepl("germline", OriginSimple, ignore.case = TRUE) |
         grepl("germline", Origin, ignore.case = TRUE) |
         is.na(OriginSimple))
  ]
  x <- x[!grepl("conflicting|association|risk factor|drug response|protective|somatic", ClinicalSignificance, ignore.case = TRUE)]
  x
}

collapse_variants <- function(x) {
  data.table::setorder(x, variant_id, -Assembly, GeneSymbol, Start)
  x[, .(
    allele_id = as.character(AlleleID[1]),
    release_date = release_date[1],
    clinical_significance = ClinicalSignificance[1],
    clin_sig_simple = ClinSigSimple[1],
    last_evaluated = LastEvaluatedDate[1],
    origin = Origin[1],
    origin_simple = OriginSimple[1],
    assembly = Assembly[1],
    variant_type = Type[1],
    variant_name = Name[1],
    gene_id = GeneID[1],
    gene_symbol = GeneSymbol[1],
    chromosome = as.character(Chromosome[1]),
    position = suppressWarnings(as.numeric(Start[1])),
    stop = suppressWarnings(as.numeric(Stop[1])),
    reference_allele = ReferenceAllele[1],
    alternative_allele = AlternateAllele[1],
    rsid = `RS# (dbSNP)`[1],
    review_status = ReviewStatus[1],
    number_submitters = suppressWarnings(as.integer(NumberSubmitters[1])),
    number_conditions = data.table::uniqueN(Name, na.rm = TRUE)
  ), by = variant_id]
}

metadata <- read.csv("data/raw/clinvar/release_metadata.csv")
baseline_meta <- metadata[metadata$release_role == "baseline_t0", ]
followup_meta <- metadata[metadata$release_role != "baseline_t0", ]

t0_raw <- read_release_dt(baseline_meta$local_file[1], as.Date(baseline_meta$release_date[1]))
t0 <- collapse_variants(t0_raw)
t0_vus_ids <- t0[clinical_significance == "Uncertain significance", variant_id]

followups <- list()
for (i in seq_len(nrow(followup_meta))) {
  release_date <- as.Date(followup_meta$release_date[i])
  release_month <- sub(".*variant_summary_([0-9]{4}-[0-9]{2}).*", "\\1", followup_meta$local_file[i])
  message("Parsing follow-up ClinVar release ", release_month)
  tx <- read_release_dt(followup_meta$local_file[i], release_date)
  tx <- tx[
    variant_id %in% t0_vus_ids &
      ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic", "Benign", "Likely benign")
  ]
  if (nrow(tx)) {
    cx <- collapse_variants(tx)
    cx$release_month <- release_month
    followups[[release_month]] <- cx
  }
}
t1_all <- data.table::rbindlist(followups, fill = TRUE)
data.table::setorder(t1_all, variant_id, release_date)
t1 <- t1_all[, .SD[1], by = variant_id]

if (nrow(t1) == 0) {
  stop("No real ClinVar VUS-to-resolved transitions found between selected releases.", call. = FALSE)
}

saveRDS(as.data.frame(t0), "data/interim/clinvar_t0_real.rds")
saveRDS(as.data.frame(t1_all), "data/interim/clinvar_followup_resolutions_all_real.rds")
saveRDS(as.data.frame(t1), "data/interim/clinvar_t1_real.rds")
utils::write.csv(head(as.data.frame(t0), 1000), "data/interim/clinvar_t0_preview.csv", row.names = FALSE)
utils::write.csv(head(as.data.frame(t1), 1000), "data/interim/clinvar_t1_preview.csv", row.names = FALSE)

log_message("outputs/logs/data_sources.txt", sprintf("Parsed real ClinVar baseline rows: %d raw, %d collapsed, %d baseline VUS.", nrow(t0_raw), nrow(t0), length(t0_vus_ids)))
log_message("outputs/logs/data_sources.txt", sprintf("Parsed real ClinVar follow-up releases and identified %d earliest VUS resolutions across %d follow-up releases.", nrow(t1), length(followups)))
