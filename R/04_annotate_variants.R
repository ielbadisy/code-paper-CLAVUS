source("R/00_setup.R")
load_or_stop("httr")
load_or_stop("jsonlite")

dat <- readRDS("data/interim/longitudinal_vus_reclassifications.rds")

annotation_cache <- "data/raw/external/myvariant_annotations.csv"

make_hgvs <- function(df) {
  ok <- df$variant_type == "single nucleotide variant" &
    !is.na(df$chromosome) & !is.na(df$position) &
    !is.na(df$reference_allele) & !is.na(df$alternative_allele) &
    nchar(df$reference_allele) == 1 & nchar(df$alternative_allele) == 1
  out <- rep(NA_character_, nrow(df))
  out[ok] <- paste0("chr", df$chromosome[ok], ":g.", as.integer(df$position[ok]),
                    df$reference_allele[ok], ">", df$alternative_allele[ok])
  out
}

pick_num <- function(x, path) {
  cur <- x
  for (p in path) {
    if (is.null(cur[[p]])) return(NA_real_)
    cur <- cur[[p]]
  }
  if (is.list(cur)) cur <- unlist(cur, use.names = FALSE)
  suppressWarnings(as.numeric(cur[1]))
}

empty_annotation_row <- function(query) {
  data.frame(
    hgvs_query = query, myvariant_found = FALSE,
    gnomad_genome_af = NA_real_, gnomad_exome_af = NA_real_,
    gnomad_genome_ac = NA_real_, gnomad_exome_ac = NA_real_,
    gnomad_genome_an = NA_real_, gnomad_exome_an = NA_real_,
    gnomad_genome_hom = NA_real_, gnomad_exome_hom = NA_real_,
    cadd_score = NA_real_, revel_score = NA_real_, sift_score = NA_real_,
    polyphen_score = NA_real_, mutationtaster_score = NA_real_,
    alphamissense_score = NA_real_, spliceai_score = NA_real_,
    stringsAsFactors = FALSE
  )
}

query_myvariant <- function(queries, batch_size = 100) {
  fields <- paste(c(
    "gnomad_genome.af.af", "gnomad_exome.af.af",
    "gnomad_genome.ac.ac", "gnomad_exome.ac.ac",
    "gnomad_genome.an.an", "gnomad_exome.an.an",
    "gnomad_genome.hom.hom", "gnomad_exome.hom.hom",
    "cadd.phred", "dbnsfp.revel.score", "dbnsfp.sift.score",
    "dbnsfp.polyphen2.hdiv.score", "dbnsfp.mutationtaster.score",
    "dbnsfp.alphamissense.score", "dbnsfp.spliceai.ds_ag",
    "dbnsfp.spliceai.ds_al", "dbnsfp.spliceai.ds_dg", "dbnsfp.spliceai.ds_dl"
  ), collapse = ",")
  chunks <- split(queries, ceiling(seq_along(queries) / batch_size))
  rows <- list()
  for (i in seq_along(chunks)) {
    message(sprintf("MyVariant batch %d/%d", i, length(chunks)))
    resp <- NULL
    for (attempt in 1:4) {
      resp <- tryCatch(
        httr::POST(
          "https://myvariant.info/v1/query",
          body = list(q = paste(chunks[[i]], collapse = ","), scopes = "_id", fields = fields),
          encode = "form",
          httr::timeout(120)
        ),
        error = function(e) e
      )
      if (!inherits(resp, "error") && httr::status_code(resp) == 200) break
      Sys.sleep(attempt)
    }
    if (inherits(resp, "error") || httr::status_code(resp) != 200) {
      warning(sprintf("MyVariant batch %d failed; marking %d variants as externally unannotated.", i, length(chunks[[i]])))
      rows <- c(rows, lapply(chunks[[i]], empty_annotation_row))
      next
    }
    parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
    rows <- c(rows, lapply(parsed, function(z) {
      splice_vals <- c(
        pick_num(z, c("dbnsfp", "spliceai", "ds_ag")),
        pick_num(z, c("dbnsfp", "spliceai", "ds_al")),
        pick_num(z, c("dbnsfp", "spliceai", "ds_dg")),
        pick_num(z, c("dbnsfp", "spliceai", "ds_dl"))
      )
      data.frame(
        hgvs_query = z$query %||% NA_character_,
        myvariant_found = is.null(z$notfound) || !isTRUE(z$notfound),
        gnomad_genome_af = pick_num(z, c("gnomad_genome", "af", "af")),
        gnomad_exome_af = pick_num(z, c("gnomad_exome", "af", "af")),
        gnomad_genome_ac = pick_num(z, c("gnomad_genome", "ac", "ac")),
        gnomad_exome_ac = pick_num(z, c("gnomad_exome", "ac", "ac")),
        gnomad_genome_an = pick_num(z, c("gnomad_genome", "an", "an")),
        gnomad_exome_an = pick_num(z, c("gnomad_exome", "an", "an")),
        gnomad_genome_hom = pick_num(z, c("gnomad_genome", "hom", "hom")),
        gnomad_exome_hom = pick_num(z, c("gnomad_exome", "hom", "hom")),
        cadd_score = pick_num(z, c("cadd", "phred")),
        revel_score = pick_num(z, c("dbnsfp", "revel", "score")),
        sift_score = pick_num(z, c("dbnsfp", "sift", "score")),
        polyphen_score = pick_num(z, c("dbnsfp", "polyphen2", "hdiv", "score")),
        mutationtaster_score = pick_num(z, c("dbnsfp", "mutationtaster", "score")),
        alphamissense_score = pick_num(z, c("dbnsfp", "alphamissense", "score")),
        spliceai_score = suppressWarnings(max(splice_vals, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
    Sys.sleep(0.2)
  }
  out <- do.call(rbind, rows)
  out$spliceai_score[is.infinite(out$spliceai_score)] <- NA_real_
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

dat$hgvs_query <- make_hgvs(dat)
queries <- sort(unique(na.omit(dat$hgvs_query)))

if (!file.exists(annotation_cache)) {
  if (!length(queries)) stop("No SNV coordinate queries available for real external annotation.", call. = FALSE)
  ann <- query_myvariant(queries)
  utils::write.csv(ann, annotation_cache, row.names = FALSE)
} else {
  ann <- read.csv(annotation_cache)
}

dat <- merge(dat, ann, by.x = "hgvs_query", by.y = "hgvs_query", all.x = TRUE, sort = FALSE)
dat$gnomad_af <- pmax(dat$gnomad_genome_af, dat$gnomad_exome_af, na.rm = TRUE)
dat$gnomad_af[is.infinite(dat$gnomad_af)] <- NA_real_
dat$gnomad_ac <- pmax(dat$gnomad_genome_ac, dat$gnomad_exome_ac, na.rm = TRUE)
dat$gnomad_ac[is.infinite(dat$gnomad_ac)] <- NA_real_
dat$gnomad_an <- pmax(dat$gnomad_genome_an, dat$gnomad_exome_an, na.rm = TRUE)
dat$gnomad_an[is.infinite(dat$gnomad_an)] <- NA_real_
dat$gnomad_hom_count <- pmax(dat$gnomad_genome_hom, dat$gnomad_exome_hom, na.rm = TRUE)
dat$gnomad_hom_count[is.infinite(dat$gnomad_hom_count)] <- NA_real_
dat$dbnsfp_available <- rowSums(!is.na(dat[, c("cadd_score", "revel_score", "sift_score", "polyphen_score",
                                               "mutationtaster_score", "alphamissense_score", "spliceai_score")])) > 0

coverage <- data.frame(
  annotation = c("myvariant_found", "gnomad_af", "cadd_score", "revel_score", "alphamissense_score", "spliceai_score", "any_dbnsfp_score"),
  n_available = c(
    sum(dat$myvariant_found %in% TRUE, na.rm = TRUE),
    sum(!is.na(dat$gnomad_af)),
    sum(!is.na(dat$cadd_score)),
    sum(!is.na(dat$revel_score)),
    sum(!is.na(dat$alphamissense_score)),
    sum(!is.na(dat$spliceai_score)),
    sum(dat$dbnsfp_available, na.rm = TRUE)
  ),
  n_total = nrow(dat)
)
coverage$pct_available <- round(100 * coverage$n_available / coverage$n_total, 2)
utils::write.csv(coverage, "outputs/tables/external_annotation_coverage.csv", row.names = FALSE)

warnings <- c(
  "Real external annotations are queried from MyVariant.info using GRCh38-style coordinate HGVS for SNVs.",
  "External annotation coverage is incomplete and is quantified in outputs/tables/external_annotation_coverage.csv.",
  "No external score is simulated or backfilled when unavailable."
)
writeLines(warnings, "outputs/logs/annotation_warnings.txt")

saveRDS(dat, "data/interim/annotated_vus_reclassifications.rds")
