# ============================================================
# RADIOMAP-IvyGAP: Data Download and Audit Pipeline
# Script 00: Download both datasets, match patient IDs, GO/NO-GO gate
# ============================================================
# Author: Daniele Piccolo, MD
# Date: 2026-02-21
# Dependencies: BiocManager, ivygapSE, SummarizedExperiment,
#               dplyr, tidyr, readxl
# Output:
#   data/raw/rnaseq/ivygap_expression_matrix.rds
#   data/raw/rnaseq/ivygap_sample_metadata.rds
#   data/raw/rnaseq/ivygap_tumor_details.rds
#   data/raw/radiomics/radiomics_features_raw.rds
#   data/processed/matched_patients.csv
#   data/processed/zone_sample_counts.csv
#   data/processed/data_audit_report.txt
# ============================================================

cat("============================================================\n")
cat("RADIOMAP-IvyGAP: Script 00 — Download and Audit\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("============================================================\n\n")

# ============================================================
# SECTION 1: SETUP
# ============================================================

cat("--- Section 1: Setup ---\n")

PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RAW_DIR     <- file.path(DATA_DIR, "raw")
PROC_DIR    <- file.path(DATA_DIR, "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RAW_RNASEQ  <- file.path(RAW_DIR, "rnaseq")
RAW_RADIO   <- file.path(RAW_DIR, "radiomics")

for (d in c(DATA_DIR, RAW_DIR, PROC_DIR, FIG_DIR, RAW_RNASEQ, RAW_RADIO)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
cat("  Directory structure created.\n")

# Install/load packages (conditional install)
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

install_if_missing("BiocManager")
install_if_missing("ivygapSE", bioc = TRUE)
install_if_missing("SummarizedExperiment", bioc = TRUE)
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("readxl")

cat("  All packages loaded.\n\n")

# ============================================================
# SECTION 2: DOWNLOAD IvyGAP RNA-seq
# ============================================================

cat("--- Section 2: Download IvyGAP RNA-seq ---\n")

expr_file  <- file.path(RAW_RNASEQ, "ivygap_expression_matrix.rds")
meta_file  <- file.path(RAW_RNASEQ, "ivygap_sample_metadata.rds")
tumor_file <- file.path(RAW_RNASEQ, "ivygap_tumor_details.rds")

if (file.exists(expr_file) && file.exists(meta_file)) {
  cat("  RNA-seq data already exists. Loading from cache.\n")
  expr_matrix <- readRDS(expr_file)
  sample_meta <- readRDS(meta_file)
  if (file.exists(tumor_file)) tumor_details <- readRDS(tumor_file)
} else {

  # --- Primary: ivygapSE Bioconductor package ---
  cat("  Loading ivygapSE Bioconductor package...\n")
  ivySE_loaded <- tryCatch({
    data(ivySE, package = "ivygapSE")
    TRUE
  }, error = function(e) {
    cat(sprintf("  WARNING: ivygapSE load failed: %s\n", e$message))
    FALSE
  })

  if (ivySE_loaded) {
    # Extract FPKM expression matrix
    expr_matrix <- SummarizedExperiment::assay(ivySE)
    cat(sprintf("  Expression matrix: %d genes x %d samples\n",
                nrow(expr_matrix), ncol(expr_matrix)))

    # Extract sample metadata
    sample_meta <- as.data.frame(SummarizedExperiment::colData(ivySE))
    cat(sprintf("  Sample metadata: %d samples, %d columns\n",
                nrow(sample_meta), ncol(sample_meta)))
    cat(sprintf("  Columns: %s\n",
                paste(colnames(sample_meta), collapse = ", ")))

    # Extract tumor/clinical details
    # metadata() is from S4Vectors, not SummarizedExperiment
    md <- S4Vectors::metadata(ivySE)
    if ("tumorDetails" %in% names(md)) {
      tumor_details <- md$tumorDetails
      cat(sprintf("  Tumor details: %d entries, columns: %s\n",
                  nrow(tumor_details),
                  paste(colnames(tumor_details), collapse = ", ")))
      saveRDS(tumor_details, tumor_file)
    } else {
      cat("  No tumorDetails in metadata. Listing available metadata:\n")
      cat(sprintf("    %s\n", paste(names(md), collapse = ", ")))
    }

    # Save
    saveRDS(expr_matrix, expr_file)
    saveRDS(sample_meta, meta_file)
    cat("  RNA-seq data saved to data/raw/rnaseq/\n")

  } else {
    # --- Fallback: Direct HTTP from Allen Institute ---
    cat("  Fallback: downloading directly from Allen Institute API...\n")

    fpkm_url  <- "http://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915"
    meta_url  <- "https://glioblastoma.alleninstitute.org/api/v2/gbm/rna_seq_samples_details.csv"
    tumor_url <- "https://glioblastoma.alleninstitute.org/api/v2/gbm/tumor_details.csv"

    # Download FPKM ZIP
    fpkm_zip <- file.path(RAW_RNASEQ, "fpkm_table.zip")
    tryCatch({
      download.file(fpkm_url, fpkm_zip, mode = "wb", quiet = FALSE)
      cat("  FPKM ZIP downloaded.\n")
    }, error = function(e) {
      stop("Failed to download FPKM data from Allen Institute.\n",
           "  Please install ivygapSE manually: BiocManager::install('ivygapSE')\n",
           "  Error: ", e$message)
    })

    # Unzip and read FPKM
    fpkm_dir <- file.path(RAW_RNASEQ, "fpkm_unzipped")
    dir.create(fpkm_dir, showWarnings = FALSE)
    unzip(fpkm_zip, exdir = fpkm_dir)
    fpkm_files <- list.files(fpkm_dir, pattern = "\\.csv$", full.names = TRUE,
                             recursive = TRUE)
    if (length(fpkm_files) == 0) {
      # Might be a TSV or TXT
      fpkm_files <- list.files(fpkm_dir, full.names = TRUE, recursive = TRUE)
    }
    cat(sprintf("  FPKM files found: %s\n", paste(basename(fpkm_files), collapse = ", ")))
    expr_df <- read.csv(fpkm_files[1], stringsAsFactors = FALSE, check.names = FALSE)
    # Assume first column is gene identifier, rest are sample columns
    gene_ids <- expr_df[[1]]
    expr_matrix <- as.matrix(expr_df[, -1])
    rownames(expr_matrix) <- gene_ids
    cat(sprintf("  Expression matrix: %d genes x %d samples\n",
                nrow(expr_matrix), ncol(expr_matrix)))

    # Download sample metadata
    meta_dest <- file.path(RAW_RNASEQ, "rna_seq_samples_details.csv")
    download.file(meta_url, meta_dest, mode = "wb", quiet = TRUE)
    sample_meta <- read.csv(meta_dest, stringsAsFactors = FALSE)
    cat(sprintf("  Sample metadata: %d rows, columns: %s\n",
                nrow(sample_meta),
                paste(colnames(sample_meta), collapse = ", ")))

    # Download tumor details
    tumor_dest <- file.path(RAW_RNASEQ, "tumor_details.csv")
    tryCatch({
      download.file(tumor_url, tumor_dest, mode = "wb", quiet = TRUE)
      tumor_details <- read.csv(tumor_dest, stringsAsFactors = FALSE)
      saveRDS(tumor_details, tumor_file)
    }, error = function(e) {
      cat(sprintf("  WARNING: Could not download tumor details: %s\n", e$message))
    })

    # Save parsed data
    saveRDS(expr_matrix, expr_file)
    saveRDS(sample_meta, meta_file)
    cat("  Fallback data saved to data/raw/rnaseq/\n")
  }
}

# Validate RNA-seq data
stopifnot(
  "Expression matrix too small (< 20,000 genes)" = nrow(expr_matrix) >= 20000,
  "Too few samples (< 200)" = ncol(expr_matrix) >= 200
)
cat(sprintf("  VALIDATED: %d genes, %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# Inspect zone/structure column
zone_col <- NULL
for (candidate in c("structure_acronym", "structure_abbreviation",
                     "structure_name", "structure_color")) {
  if (candidate %in% colnames(sample_meta)) {
    zone_col <- candidate
    break
  }
}
if (is.null(zone_col)) {
  cat("  WARNING: No zone column found. Available columns:\n")
  cat(sprintf("    %s\n", paste(colnames(sample_meta), collapse = ", ")))
  # Try to find any column with CT/IT/LE values
  for (col in colnames(sample_meta)) {
    vals <- unique(as.character(sample_meta[[col]]))
    if (any(grepl("^CT$|^IT$|^LE$|Cellular Tumor|Leading Edge", vals))) {
      zone_col <- col
      cat(sprintf("  Found zone column by content: %s\n", col))
      break
    }
  }
}
cat(sprintf("  Zone column: %s\n", zone_col))
cat(sprintf("  Zone values: %s\n",
            paste(sort(unique(as.character(sample_meta[[zone_col]]))), collapse = ", ")))

# Inspect tumor/patient ID column
tumor_col <- NULL
for (candidate in c("tumor_name", "tumor_id", "donor_name", "donor_id",
                     "specimen_name")) {
  if (candidate %in% colnames(sample_meta)) {
    tumor_col <- candidate
    break
  }
}
cat(sprintf("  Tumor ID column: %s\n", tumor_col))
cat(sprintf("  First 5 tumor IDs: %s\n",
            paste(head(unique(as.character(sample_meta[[tumor_col]])), 5),
                  collapse = ", ")))

# Extract base zone from compound structure_acronym labels
# e.g., "CT-CD44" -> "CT", "CTmvp-ITGA6" -> "CTmvp", "IT-reference-histology" -> "IT"
# The base zone is always the first token when splitting on "-"
sample_meta$zone <- sapply(
  strsplit(as.character(sample_meta[[zone_col]]), "-"),
  `[`, 1
)
cat(sprintf("  Base zones extracted: %s\n",
            paste(sort(unique(sample_meta$zone)), collapse = ", ")))

# Extract patient W-number from tumor identifier
# Pattern: "W31" from entries like "W31", "W31-1-1", "W31-1-1-E.1.3"
sample_meta$patient_id <- sub("^(W\\d+).*", "\\1",
                              as.character(sample_meta[[tumor_col]]))
n_patients_rnaseq <- length(unique(sample_meta$patient_id))
cat(sprintf("  Unique patients (RNA-seq): %d\n", n_patients_rnaseq))
cat(sprintf("  Patient IDs: %s\n",
            paste(sort(unique(sample_meta$patient_id)), collapse = ", ")))
cat("\n")

# ============================================================
# SECTION 3: DOWNLOAD IVYGAP-RADIOMICS
# ============================================================

cat("--- Section 3: Download IVYGAP-RADIOMICS ---\n")

radio_rds <- file.path(RAW_RADIO, "radiomics_features_raw.rds")

if (file.exists(radio_rds)) {
  cat("  Radiomics data already parsed. Loading from cache.\n")
  radio_all <- readRDS(radio_rds)
} else {

  # Check if ZIP or extracted files already exist
  zip_files <- list.files(RAW_RADIO, pattern = "\\.zip$", full.names = TRUE)

  if (length(zip_files) == 0) {
    cat("  Attempting download from TCIA wiki...\n")

    # TCIA wiki attachment URLs (page ID = 70222827)
    # Actual filenames discovered via Confluence REST API
    tcia_urls <- c(
      "https://wiki.cancerimagingarchive.net/download/attachments/70222827/Multi-Institutional%20Paired%20Expert%20Segmentations%20Radiomic%20Features%20%20and%20Reproducibility%20Evaluation%20on%20SRI.zip?api=v2",
      "https://wiki.cancerimagingarchive.net/download/attachments/70222827/Multi-Institutional%20Paired%20Expert%20Segmentations%20Radiomic%20Features%20and%20Reproducibility%20Evaluation%20on%20SRI.zip?api=v2"
    )

    zip_dest <- file.path(RAW_RADIO, "ivygap_radiomics.zip")
    downloaded <- FALSE

    for (url in tcia_urls) {
      tryCatch({
        download.file(url, zip_dest, mode = "wb", quiet = TRUE)
        # Verify it's a real ZIP (> 10 KB, starts with PK signature)
        if (file.exists(zip_dest) && file.size(zip_dest) > 10000) {
          con <- file(zip_dest, "rb")
          sig <- readBin(con, "raw", 2)
          close(con)
          if (identical(sig, charToRaw("PK"))) {
            downloaded <- TRUE
            cat(sprintf("  Downloaded from: %s\n", url))
            break
          }
        }
        if (file.exists(zip_dest) && !downloaded) file.remove(zip_dest)
      }, error = function(e) {
        if (file.exists(zip_dest)) file.remove(zip_dest)
      })
    }

    if (!downloaded) {
      # Try TCIA metadata CSV separately
      meta_csv_url <- "https://wiki.cancerimagingarchive.net/download/attachments/70222827/ivygap_metadata.csv?api=v2"
      tryCatch({
        download.file(meta_csv_url,
                      file.path(RAW_RADIO, "tcia_metadata.csv"),
                      mode = "wb", quiet = TRUE)
        cat("  Downloaded TCIA metadata CSV.\n")
      }, error = function(e) NULL)

      stop(paste0(
        "\n=== MANUAL DOWNLOAD REQUIRED ===\n",
        "Automatic download of IVYGAP-RADIOMICS failed.\n",
        "Please download manually:\n",
        "  1. Go to: https://www.cancerimagingarchive.net/analysis-result/ivygap-radiomics/\n",
        "  2. Download the 'Radiomic Features & Reproducibility' ZIP (~17 MB)\n",
        "  3. Place the ZIP file in:\n     ", RAW_RADIO, "\n",
        "  4. Re-run this script.\n"
      ))
    }

    zip_files <- zip_dest
  }

  # Unzip
  cat("  Unzipping radiomics data...\n")
  unzip_dir <- file.path(RAW_RADIO, "unzipped")
  dir.create(unzip_dir, showWarnings = FALSE, recursive = TRUE)
  unzip(zip_files[1], exdir = unzip_dir)

  # List contents
  all_files <- list.files(unzip_dir, recursive = TRUE, full.names = TRUE)
  cat(sprintf("  ZIP contents (%d files):\n", length(all_files)))
  for (f in all_files) {
    cat(sprintf("    %s (%.1f KB)\n", basename(f), file.size(f) / 1024))
  }

  # Find and read all tabular files (CSV, XLSX)
  csv_files  <- grep("\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
  xlsx_files <- grep("\\.xlsx?$", all_files, value = TRUE, ignore.case = TRUE)
  tabular_files <- c(csv_files, xlsx_files)

  if (length(tabular_files) == 0) {
    stop("No CSV or XLSX files found in the radiomics ZIP.\n",
         "  Files found: ", paste(basename(all_files), collapse = ", "))
  }

  # Read all tabular files, store in a named list
  radio_all <- list()
  for (f in tabular_files) {
    nm <- tools::file_path_sans_ext(basename(f))
    tryCatch({
      if (grepl("\\.csv$", f, ignore.case = TRUE)) {
        df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        # For XLSX, read all sheets
        sheets <- readxl::excel_sheets(f)
        if (length(sheets) == 1) {
          df <- as.data.frame(readxl::read_excel(f))
        } else {
          # Read each sheet as a separate entry
          for (s in sheets) {
            sdf <- as.data.frame(readxl::read_excel(f, sheet = s))
            radio_all[[paste0(nm, "_", gsub(" ", "_", s))]] <- sdf
            cat(sprintf("    Read: %s [sheet: %s] (%d x %d)\n",
                        basename(f), s, nrow(sdf), ncol(sdf)))
          }
          next
        }
      }
      radio_all[[nm]] <- df
      cat(sprintf("    Read: %s (%d x %d)\n", basename(f), nrow(df), ncol(df)))
    }, error = function(e) {
      cat(sprintf("    Could not read %s: %s\n", basename(f), e$message))
    })
  }

  # Save all parsed radiomics data
  saveRDS(radio_all, radio_rds)
  cat("  Radiomics data saved to data/raw/radiomics/\n")
}

# --- Restructure radiomics from Combined.xlsx sheets ---
# Data structure: 12 sheets in Combined.xlsx, one per {sequence}_{subcompartment}
# Each sheet: 2360 features (rows) x 32 cols (1 "Features" col + 31 patient W-numbers)
# Subcompartments: et (Enhancing Tumor), nec (Necrosis/NET), ed (Edema)
# Sequences: T1, T2, T1ce, FLAIR

cat("\n  Restructuring radiomics from Combined sheets...\n")

# Map sheet suffixes to standard subcompartment names
subcomp_map <- c(et = "ET", nec = "NET", ed = "ED")

# Find all Combined_*_* sheets (pattern: Combined_{seq}_{subcomp})
combined_keys <- grep("^Combined_", names(radio_all), value = TRUE)
cat(sprintf("  Found %d Combined sheets: %s\n", length(combined_keys),
            paste(combined_keys, collapse = ", ")))

if (length(combined_keys) == 0) {
  stop("No Combined_* sheets found in radiomics data. ",
       "Expected sheets like Combined_T1_et, Combined_T2_nec, etc.")
}

# Build per-subcompartment feature matrices (patients as rows, features as columns)
radio_per_sc <- list()

for (sc_short in names(subcomp_map)) {
  sc_label <- subcomp_map[sc_short]
  sc_sheets <- grep(paste0("_", sc_short, "$"), combined_keys, value = TRUE)

  if (length(sc_sheets) == 0) {
    cat(sprintf("  WARNING: No sheets found for subcompartment %s\n", sc_label))
    next
  }

  cat(sprintf("  %s (%s): combining %d sequence sheets\n",
              sc_label, sc_short, length(sc_sheets)))

  # Combine features across MRI sequences
  feat_list <- list()
  for (sheet_key in sc_sheets) {
    df <- radio_all[[sheet_key]]
    # Column 1 = "Features" (feature names), rest = patient W-number columns
    feat_names <- as.character(df[[1]])
    patient_cols <- colnames(df)[-1]

    # Extract sequence name from sheet key (Combined_{seq}_{subcomp})
    seq_name <- sub("^Combined_(.+)_[^_]+$", "\\1", sheet_key)

    # Prefix feature names with sequence to ensure uniqueness
    prefixed_names <- paste0(seq_name, "_", feat_names)

    # Transpose: patients become rows, features become columns
    feat_mat <- t(as.matrix(df[, -1, drop = FALSE]))
    colnames(feat_mat) <- prefixed_names
    # Convert to numeric (in case of character issues from Excel)
    storage.mode(feat_mat) <- "double"

    feat_list[[seq_name]] <- feat_mat
    cat(sprintf("    %s: %d features\n", seq_name, ncol(feat_mat)))
  }

  # cbind across sequences (same patients, different features)
  combined_mat <- do.call(cbind, feat_list)
  sc_df <- as.data.frame(combined_mat, stringsAsFactors = FALSE)
  sc_df$patient_id <- rownames(combined_mat)
  radio_per_sc[[sc_label]] <- sc_df

  cat(sprintf("    %s total: %d patients x %d features\n",
              sc_label, nrow(sc_df), ncol(sc_df) - 1))
}

# Extract patient IDs from any subcompartment
radio_patient_ids <- sort(unique(radio_per_sc[[1]]$patient_id))
n_subjects_radio <- length(radio_patient_ids)

cat(sprintf("  Unique radiomics subjects: %d\n", n_subjects_radio))
cat(sprintf("  Radiomics patient IDs: %s\n",
            paste(radio_patient_ids, collapse = ", ")))

# Save structured radiomics data for script 02
saveRDS(radio_per_sc, file.path(PROC_DIR, "radiomics_per_subcompartment.rds"))
cat("  Structured radiomics saved to data/processed/\n")

data_format <- "sheets"
id_col_radio <- "patient_id"
cat(sprintf("  Data format: %s (12 sheets restructured)\n\n", data_format))

# ============================================================
# SECTION 4: PATIENT ID MATCHING
# ============================================================

cat("--- Section 4: Patient ID Matching ---\n")

rnaseq_patients <- sort(unique(sample_meta$patient_id))
radio_patients  <- sort(unique(radio_patient_ids))

matched_patients <- intersect(rnaseq_patients, radio_patients)
matched_N <- length(matched_patients)

cat(sprintf("  RNA-seq patients: %d\n", length(rnaseq_patients)))
cat(sprintf("  Radiomics patients: %d\n", length(radio_patients)))
cat(sprintf("  MATCHED patients: %d\n", matched_N))

if (matched_N > 0) {
  cat(sprintf("  Matched IDs: %s\n", paste(sort(matched_patients), collapse = ", ")))
}

# Show unmatched for debugging
rnaseq_only <- setdiff(rnaseq_patients, radio_patients)
radio_only  <- setdiff(radio_patients, rnaseq_patients)
if (length(rnaseq_only) > 0) {
  cat(sprintf("  RNA-seq only (%d): %s\n", length(rnaseq_only),
              paste(head(rnaseq_only, 10), collapse = ", ")))
}
if (length(radio_only) > 0) {
  cat(sprintf("  Radiomics only (%d): %s\n", length(radio_only),
              paste(head(radio_only, 10), collapse = ", ")))
}
cat("\n")

# ============================================================
# SECTION 5: GO/NO-GO GATE
# ============================================================

cat("--- Section 5: GO/NO-GO Gate ---\n")
cat(sprintf("  Matched N = %d (threshold = 25)\n", matched_N))

if (matched_N < 25) {
  stop(sprintf(
    paste0("\n",
           "========================================\n",
           "  KILL CRITERION MET: matched N = %d < 25\n",
           "  PROJECT TERMINATED\n",
           "========================================\n"),
    matched_N
  ))
}

cat(sprintf("  *** GO: %d matched patients (>= 25 threshold) ***\n\n", matched_N))

# ============================================================
# SECTION 6: ZONE AVAILABILITY AUDIT
# ============================================================

cat("--- Section 6: Zone Availability Audit ---\n")

# Subset to matched patients
matched_meta <- sample_meta[sample_meta$patient_id %in% matched_patients, ]

# Define zones of interest (for primary mapping)
primary_zones <- c("CT", "CTmvp", "CTpan", "IT", "LE")

# Tabulate samples per zone per patient
zone_counts <- matched_meta %>%
  filter(zone %in% primary_zones) %>%
  group_by(patient_id, zone) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = zone,
    values_from = n_samples,
    values_fill = 0
  )

cat(sprintf("  Matched patients with zone data: %d\n", nrow(zone_counts)))
cat("  Samples per zone per patient:\n")

# Summary stats per zone
for (z in primary_zones) {
  if (z %in% colnames(zone_counts)) {
    vals <- zone_counts[[z]]
    n_with <- sum(vals > 0)
    cat(sprintf("    %s: %d/%d patients have data (median %.0f, range %d-%d samples)\n",
                z, n_with, nrow(zone_counts), median(vals[vals > 0]),
                min(vals[vals > 0]), max(vals[vals > 0])))
  } else {
    cat(sprintf("    %s: NO DATA\n", z))
  }
}

# Check for complete zone coverage (all 5 primary zones)
if (all(primary_zones %in% colnames(zone_counts))) {
  complete_coverage <- zone_counts %>%
    filter(CT > 0, CTmvp > 0, CTpan > 0, IT > 0, LE > 0)
  cat(sprintf("  Patients with ALL 5 primary zones: %d/%d\n",
              nrow(complete_coverage), nrow(zone_counts)))
} else {
  cat("  WARNING: Not all primary zones found in data.\n")
}

# Check subcompartment coverage (what matters for analysis)
# ET needs CT or CTmvp; NET needs CTpan; ED needs IT or LE
subcomp_coverage <- zone_counts %>%
  mutate(
    has_ET  = (CT > 0 | CTmvp > 0),
    has_NET = (CTpan > 0),
    has_ED  = (IT > 0 | LE > 0),
    n_subcomp = has_ET + has_NET + has_ED
  )
cat(sprintf("  Subcompartment coverage:\n"))
cat(sprintf("    ET (CT|CTmvp): %d/%d patients\n",
            sum(subcomp_coverage$has_ET), nrow(subcomp_coverage)))
cat(sprintf("    NET (CTpan): %d/%d patients\n",
            sum(subcomp_coverage$has_NET), nrow(subcomp_coverage)))
cat(sprintf("    ED (IT|LE): %d/%d patients\n",
            sum(subcomp_coverage$has_ED), nrow(subcomp_coverage)))
cat(sprintf("    All 3 subcompartments: %d/%d patients\n",
            sum(subcomp_coverage$n_subcomp == 3), nrow(subcomp_coverage)))

# Also report excluded zones
all_zones <- sort(unique(matched_meta$zone))
excluded_zones <- setdiff(all_zones, primary_zones)
if (length(excluded_zones) > 0) {
  cat(sprintf("  Zones present but excluded from primary mapping: %s\n",
              paste(excluded_zones, collapse = ", ")))
  for (z in excluded_zones) {
    n_samples <- sum(matched_meta$zone == z)
    cat(sprintf("    %s: %d samples\n", z, n_samples))
  }
}
cat("\n")

# ============================================================
# SECTION 7: RADIOMIC FEATURE AUDIT
# ============================================================

cat("--- Section 7: Radiomic Feature Audit ---\n")

# Audit per-subcompartment feature matrices
feature_cols_radio <- list()
total_na <- 0
total_inf <- 0
total_cells <- 0

for (sc in names(radio_per_sc)) {
  sc_df <- radio_per_sc[[sc]]
  feat_cols <- setdiff(colnames(sc_df), "patient_id")
  feature_cols_radio[[sc]] <- feat_cols

  n_feat <- length(feat_cols)
  n_pat  <- nrow(sc_df)
  feat_mat <- as.matrix(sc_df[, feat_cols, drop = FALSE])
  n_na  <- sum(is.na(feat_mat))
  n_inf <- sum(is.infinite(feat_mat))

  total_na    <- total_na + n_na
  total_inf   <- total_inf + n_inf
  total_cells <- total_cells + length(feat_mat)

  # Count per MRI sequence
  t1_n    <- sum(grepl("^T1_", feat_cols))
  t2_n    <- sum(grepl("^T2_", feat_cols))
  t1ce_n  <- sum(grepl("^T1ce_", feat_cols))
  flair_n <- sum(grepl("^FLAIR_", feat_cols))

  cat(sprintf("  %s: %d patients x %d features\n", sc, n_pat, n_feat))
  cat(sprintf("    Per sequence: T1=%d, T2=%d, T1ce=%d, FLAIR=%d\n",
              t1_n, t2_n, t1ce_n, flair_n))
  cat(sprintf("    NA: %d (%.2f%%), Inf: %d (%.2f%%)\n",
              n_na, 100 * n_na / length(feat_mat),
              n_inf, 100 * n_inf / length(feat_mat)))
}

cat(sprintf("  Total across subcompartments: NA=%d, Inf=%d / %d cells\n",
            total_na, total_inf, total_cells))
cat("\n")

# ============================================================
# SECTION 8: SUMMARY REPORT
# ============================================================

cat("--- Section 8: Summary Report ---\n")

# Build report
report <- c(
  "============================================================",
  "RADIOMAP-IvyGAP: Data Audit Report",
  sprintf("Generated: %s", Sys.time()),
  "============================================================",
  "",
  "--- DATASETS ---",
  sprintf("IvyGAP RNA-seq (ivygapSE): %d genes x %d samples", nrow(expr_matrix), ncol(expr_matrix)),
  sprintf("IvyGAP RNA-seq patients: %d", n_patients_rnaseq),
  sprintf("IVYGAP-RADIOMICS: 12 sheets (4 sequences x 3 subcompartments), %d features/seq",
          nrow(radio_all[[combined_keys[1]]])),
  sprintf("IVYGAP-RADIOMICS subjects: %d", n_subjects_radio),
  sprintf("IVYGAP-RADIOMICS features per subcompartment: %d",
          length(feature_cols_radio[[1]])),
  sprintf("Radiomics data format: %s", data_format),
  "",
  "--- PATIENT MATCHING ---",
  sprintf("RNA-seq patients: %d", length(rnaseq_patients)),
  sprintf("Radiomics patients: %d", length(radio_patients)),
  sprintf("MATCHED: %d", matched_N),
  sprintf("Matched IDs: %s", paste(sort(matched_patients), collapse = ", ")),
  "",
  sprintf("GO/NO-GO: %s (threshold = 25)", ifelse(matched_N >= 25, "GO", "NO-GO")),
  "",
  "--- ZONE AVAILABILITY (matched patients) ---"
)

for (z in primary_zones) {
  if (z %in% colnames(zone_counts)) {
    vals <- zone_counts[[z]]
    n_with <- sum(vals > 0)
    report <- c(report,
      sprintf("  %s: %d/%d patients, median %.0f samples",
              z, n_with, nrow(zone_counts), median(vals[vals > 0])))
  }
}

report <- c(report,
  sprintf("  Patients with all 3 subcompartments: %d/%d",
          sum(subcomp_coverage$n_subcomp == 3), nrow(subcomp_coverage)),
  "",
  "--- RADIOMIC FEATURES ---",
  sprintf("  Features per subcompartment: %d", length(feature_cols_radio[[1]])),
  sprintf("  Subcompartments: %s", paste(names(radio_per_sc), collapse = ", ")),
  sprintf("  Missing values: %d (%.2f%%)", total_na,
          ifelse(total_cells > 0, 100 * total_na / total_cells, 0)),
  sprintf("  Infinite values: %d (%.2f%%)", total_inf,
          ifelse(total_cells > 0, 100 * total_inf / total_cells, 0)),
  "",
  "--- SESSION INFO ---",
  sprintf("  R version: %s", R.version.string),
  sprintf("  Platform: %s", R.version$platform),
  ""
)

# Write report
report_file <- file.path(PROC_DIR, "data_audit_report.txt")
writeLines(report, report_file)
cat(sprintf("  Audit report saved to: %s\n", report_file))

# Save matched patients
matched_df <- data.frame(
  patient_id = sort(matched_patients),
  in_rnaseq  = TRUE,
  in_radiomics = TRUE,
  stringsAsFactors = FALSE
)
matched_csv <- file.path(PROC_DIR, "matched_patients.csv")
write.csv(matched_df, matched_csv, row.names = FALSE)
cat(sprintf("  Matched patients saved to: %s\n", matched_csv))

# Save zone sample counts
zone_csv <- file.path(PROC_DIR, "zone_sample_counts.csv")
write.csv(as.data.frame(zone_counts), zone_csv, row.names = FALSE)
cat(sprintf("  Zone counts saved to: %s\n", zone_csv))

# Save key variables for downstream scripts
audit_vars <- list(
  matched_patients = matched_patients,
  matched_N        = matched_N,
  zone_col         = zone_col,
  tumor_col        = tumor_col,
  id_col_radio     = id_col_radio,
  data_format      = data_format,
  feature_cols_radio = feature_cols_radio,
  primary_zones    = primary_zones
)
saveRDS(audit_vars, file.path(PROC_DIR, "audit_variables.rds"))
cat("  Audit variables saved for downstream scripts.\n")

cat("\n============================================================\n")
cat(sprintf("Script 00 complete. Matched N = %d. Status: GO.\n", matched_N))
cat("============================================================\n")
