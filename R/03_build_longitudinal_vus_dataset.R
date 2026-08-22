source("R/00_setup.R")

t0 <- readRDS("data/interim/clinvar_t0_real.rds")
t1 <- readRDS("data/interim/clinvar_t1_real.rds")

names(t0)[names(t0) == "clinical_significance"] <- "baseline_class"
names(t1)[names(t1) == "clinical_significance"] <- "final_class"
names(t0)[names(t0) == "release_date"] <- "baseline_date"
names(t1)[names(t1) == "release_date"] <- "final_date"
names(t1)[names(t1) == "review_status"] <- "final_review_status"
longitudinal <- merge(
  t0,
  t1[, c("variant_id", "final_date", "release_month", "final_class", "final_review_status")],
  by = "variant_id"
)
names(longitudinal)[names(longitudinal) == "release_month"] <- "final_release_month"

longitudinal$final_class_group <- class_group(longitudinal$final_class)
longitudinal$label <- ifelse(longitudinal$final_class_group == "P_LP", 1L,
                             ifelse(longitudinal$final_class_group == "B_LB", 0L, NA_integer_))
longitudinal$time_to_reclassification_days <- as.integer(as.Date(longitudinal$final_date) - as.Date(longitudinal$baseline_date))

included <- subset(longitudinal,
                   baseline_class == "Uncertain significance" &
                     !is.na(label) &
                     !is.na(variant_id))

if (nrow(included) < 100 || length(unique(included$label)) < 2) {
  stop(sprintf("Real ClinVar cohort is too small or single-class after VUS resolution rules: n=%d.", nrow(included)), call. = FALSE)
}

saveRDS(included, "data/interim/longitudinal_vus_reclassifications.rds")
utils::write.csv(included, "data/interim/longitudinal_vus_reclassifications.csv", row.names = FALSE)
log_message("outputs/logs/data_sources.txt", sprintf("Included %d VUS reclassifications after label rules.", nrow(included)))
