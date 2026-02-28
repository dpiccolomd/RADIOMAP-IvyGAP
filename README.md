# RADIOMAP-IvyGAP

**Subcompartment-Level Radiomic Features Associate with Regional Transcriptomic Programs in Glioblastoma: An Exploratory Analysis of the Ivy Glioblastoma Atlas Project**

Daniele Piccolo, Marco Vindigni

## Overview

This repository contains the analysis code for an exploratory study linking MRI-derived radiomic features from the IVYGAP-RADIOMICS dataset (Pati et al., 2020) with zone-specific transcriptomic pathway enrichment scores from the IvyGAP RNA-seq atlas (Puchalski et al., 2018).

The study tests whether radiomic features extracted from BraTS-style MRI subcompartments (enhancing tumor, non-enhancing tumor, peritumoral edema) associate with transcriptomic programs computed via single-sample Gene Set Enrichment Analysis (ssGSEA) using 24 gene sets across 28 matched patients.

## Repository Structure

```
RADIOMAP-IvyGAP/
├── R/
│   ├── 00_download_and_audit.R        # Data download, patient matching, GO/NO-GO gate
│   ├── 01_ssgsea_computation.R        # 24 ssGSEA pathway enrichment scores
│   ├── 02_prepare_analysis_data.R     # Zone-to-subcompartment aggregation, feature merge
│   ├── 03_statistical_analysis.R      # Mixed-effects models, Elastic Net, permutation tests
│   └── 03b_nested_cv_analysis.R       # Nested LOPO-CV (primary predictive analysis)
├── results/                           # Pre-computed output tables and model objects
├── README.md
├── .gitignore
└── LICENSE
```

## Requirements

### R (>= 4.3)

Packages are installed automatically by the scripts if missing:

- **Bioconductor**: GSVA, ivygapSE, SummarizedExperiment, msigdbr, org.Hs.eg.db, AnnotationDbi
- **CRAN**: lme4, lmerTest, performance, glmnet, pmsampsize, caret, dplyr, tidyr, broom.mixed, ggplot2, patchwork, readxl

## Reproduction

Scripts must be run **sequentially from the repository root**:

```bash
cd RADIOMAP-IvyGAP

# Step 1: Download data and match patients (~5 min, downloads ~500 MB)
Rscript R/00_download_and_audit.R

# Step 2: Compute ssGSEA enrichment scores (~10 min)
Rscript R/01_ssgsea_computation.R

# Step 3: Zone-to-subcompartment aggregation and feature merge (~1 min)
Rscript R/02_prepare_analysis_data.R

# Step 4: Mixed-effects models, Elastic Net, permutation tests (~30 min)
Rscript R/03_statistical_analysis.R

# Step 5: Nested cross-validation with internal feature selection (~3-4 hours)
#         Includes 1000-permutation test for each significant pathway
Rscript R/03b_nested_cv_analysis.R
```

**Notes:**
- Script 00 downloads both datasets automatically (IvyGAP via Bioconductor, IVYGAP-RADIOMICS from TCIA). If TCIA download fails, manual download instructions are provided.
- Script 03 saves a session file (`results/03_analysis_session.RData`, ~100 MB) that is required by script 03b.
- Script 03b runs 1000 nested permutations per significant pathway; total runtime is approximately 3-4 hours.
- The `results/` directory contains pre-computed outputs for reference. Re-running the full pipeline will regenerate these files.

## Data Sources

- **IvyGAP RNA-seq**: Allen Institute for Brain Science — [https://glioblastoma.alleninstitute.org/](https://glioblastoma.alleninstitute.org/)
- **IVYGAP-RADIOMICS**: The Cancer Imaging Archive — [https://doi.org/10.7937/9j41-7d44](https://doi.org/10.7937/9j41-7d44)

Both datasets are publicly available. No IRB approval was required.

## Key Results

| Pathway | R²_cv | Permutation p | Bootstrap 95% CI |
|---|---|---|---|
| Angiogenesis | 0.209 | 0.006 | [0.028, 0.353] |
| Inflammatory Response | 0.185 | 0.008 | [0.071, 0.355] |
| IvyGAP CTpan module | 0.133 | 0.013 | [-0.079, 0.350] |

21 of 24 pathways showed no robust radiomic signal (R²_cv ≤ 0).

## Citation

If you use this code, please cite the associated manuscript (reference to be updated upon publication).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
