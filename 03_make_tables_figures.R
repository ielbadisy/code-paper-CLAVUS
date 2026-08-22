# CLAVUS -- Step 3: Produce manuscript tables and figures
#
# Generates all LaTeX tables (booktabs/threeparttable/siunitx) and all
# publication figures (ggplot2/patchwork, ragg 400 DPI) used in the
# manuscript, including the PDP/ICE panel for the five most important
# features.
#
# Run from the project root, after 02_run_ml_pipeline.R:
#   Rscript 03_make_tables_figures.R
#
# Output: outputs/tables/tex/*.tex, outputs/figures/*.{png,pdf}

source("R/00_setup.R")
source("R/11_tables_figures.R")
source("R/11_tables_latex.R")
source("R/11_figures_production.R")
source("R/12_pdp_ice_figure.R")
