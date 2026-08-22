source("R/00_setup.R")
options(timeout = max(900, getOption("timeout")))

release_months <- data.frame(
  role = c("baseline_t0", "followup_2021", "followup_2022", "followup_2023", "followup_2024", "external_followup_2026"),
  release_month = c("2020-01", "2021-01", "2022-01", "2023-01", "2024-12", "2026-04"),
  release_date = as.Date(c("2020-01-02", "2021-01-07", "2022-01-06", "2023-01-05", "2024-12-05", "2026-04-02"))
)

archive_url <- function(month) {
  year <- substr(month, 1, 4)
  if (as.integer(year) <= 2024) {
    sprintf("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/%s/variant_summary_%s.txt.gz", year, month)
  } else {
    sprintf("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_%s.txt.gz", month)
  }
}

clinvar_releases <- transform(
  release_months,
  url = vapply(release_month, archive_url, character(1)),
  local_file = file.path("data/raw/clinvar", sprintf("variant_summary_%s.txt.gz", release_month))
)

for (i in seq_len(nrow(clinvar_releases))) {
  dest <- clinvar_releases$local_file[i]
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    log_message("outputs/logs/data_sources.txt", paste("Downloading real ClinVar release:", clinvar_releases$url[i]))
    utils::download.file(clinvar_releases$url[i], destfile = dest, mode = "wb", quiet = FALSE)
  }
}

metadata <- data.frame(
  source = "ClinVar",
  release_role = clinvar_releases$role,
  release_date = clinvar_releases$release_date,
  local_file = clinvar_releases$local_file,
  url = clinvar_releases$url,
  note = "Real ClinVar tab-delimited archived variant_summary release."
)

utils::write.csv(metadata, "data/raw/clinvar/release_metadata.csv", row.names = FALSE)
log_message("outputs/logs/data_sources.txt", "Real ClinVar archived releases prepared for multi-release reconstruction: 2020-01 baseline with 2021, 2022, 2023, 2024, and 2026 follow-up releases.", append = FALSE)
