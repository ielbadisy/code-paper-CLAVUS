# CLAVUS -- Step 1: Build the curated longitudinal VUS reclassification dataset
#
# Downloads archived ClinVar variant_summary releases, parses and links them
# across releases, queries external functional annotations (MyVariant.info),
# and engineers the final baseline predictor set.
#
# Run from the project root (this script's parent directory):
#   Rscript 01_build_dataset.R
#
# Output: data/curated/vus_resolution_dataset.rds (+ .csv) and
#         data/curated/feature_recipe.rds
# These two curated outputs are already included with this submission under
# data/curated/, so this step is optional unless re-downloading from ClinVar.

source("R/00_setup.R")
source("R/01_download_clinvar.R")
source("R/02_parse_clinvar_releases.R")
source("R/03_build_longitudinal_vus_dataset.R")
source("R/03b_download_external_features.R")
source("R/04_annotate_variants.R")
source("R/05_feature_engineering.R")
