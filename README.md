# CLAVUS -- reproducible analysis code

Complete R pipeline reproducing every number, table, and figure reported in
the manuscript "A Temporally Validated Machine Learning Framework for
Predicting Clinical Resolution of Variants of Uncertain Significance."

## Structure

```
01_build_dataset.R        entry point -- cohort construction, annotation, features
02_run_ml_pipeline.R      entry point -- learner training, selection, validation,
                           interpretability, adjusted-effect (ATE) estimation
03_make_tables_figures.R  entry point -- all manuscript tables (LaTeX) and figures

R/                        implementation scripts sourced by the three entry
                           points above, numbered in execution order (00-12)
data/raw/, data/interim/, data/curated/  empty placeholders populated by 01_build_dataset.R
outputs/                  empty placeholders populated by scripts 02 and 03
```

## Reproducing the results

The curated dataset is not redistributed in this repository. It is rebuilt
directly from the official ClinVar archive (see Data availability below):

```r
install.packages(c("funcml", "mimar"))   # both on CRAN
source("01_build_dataset.R")             # downloads ClinVar releases, builds the curated cohort
source("02_run_ml_pipeline.R")
source("03_make_tables_figures.R")
```

## Software

R 4.4.1. Key packages: `funcml` (learner training / CV / interpretability /
ATE estimation), `mimar` (chained random-forest multiple imputation),
`ggplot2`, `patchwork`, `ragg`, `knitr`, `quarto`, both available from CRAN.
Global random seed: 20260427. Exact package versions used to produce the
reported results are recorded as comments at the top of `R/00_setup.R`.

## Data availability

All data used in this study are derived from the official ClinVar archive
maintained by NCBI: `ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/`.
Raw archived `variant_summary` releases (~1.1 GB) are not redistributed here;
`01_build_dataset.R` downloads them directly from this official source and
deterministically rebuilds the curated, analysis-ready cohort
(14,046 variants) used for all reported results.

## Citation

If you use this code, please cite:

> El Badisy I, El Kadiri Y, Ghazi B, Elfahime E, Boutayeb S. CLAVUS: A
> Machine Learning Framework with Temporal Validation for Predicting
> Clinical Resolution of Variants of Uncertain Significance. (manuscript
> under review).

## License

MIT License. See [LICENSE](LICENSE).
