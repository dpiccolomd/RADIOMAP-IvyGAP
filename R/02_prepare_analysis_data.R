# ============================================================
# RADIOMAP-IvyGAP: Prepare Analysis-Ready Dataset
# Script 02: Zone-to-subcompartment aggregation + radiomics merge
# ============================================================
# Author: Daniele Piccolo, MD
# Date: 2026-02-21
# Dependencies: dplyr, tidyr
# Input:
#   data/processed/ssgsea_scores.rds (from script 01)
#   data/processed/sample_metadata_matched.rds (from script 01)
#   data/processed/audit_variables.rds (from script 00)
#   data/raw/radiomics/radiomics_features_raw.rds (from script 00)
# Output:
#   data/analysis_data.rds (input for script 03)
# ============================================================

cat("============================================================\n")
cat("RADIOMAP-IvyGAP: Script 02 — Prepare Analysis Data\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("============================================================\n\n")

# --- Setup ---
PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RAW_DIR     <- file.path(DATA_DIR, "raw")
PROC_DIR    <- file.path(DATA_DIR, "processed")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# LOAD UPSTREAM OUTPUTS
# ============================================================

cat("--- Loading upstream data ---\n")

ssgsea_scores <- readRDS(file.path(PROC_DIR, "ssgsea_scores.rds"))
meta_matched  <- readRDS(file.path(PROC_DIR, "sample_metadata_matched.rds"))
audit_vars    <- readRDS(file.path(PROC_DIR, "audit_variables.rds"))

matched_patients <- audit_vars$matched_patients

cat(sprintf("  ssGSEA scores: %d pathways x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))
cat(sprintf("  Matched patients: %d\n", length(matched_patients)))

# ============================================================
# SECTION 1: ZONE-TO-SUBCOMPARTMENT AGGREGATION (ssGSEA)
# ============================================================

cat("\n--- Section 1: Zone-to-Subcompartment Aggregation ---\n")

# Primary mapping:
#   ET  = mean(CT samples, CTmvp samples) per patient
#   NET = mean(CTpan samples) per patient
#   ED  = mean(IT samples, LE samples) per patient

zone_to_subcomp <- list(
  ET  = c("CT", "CTmvp"),
  NET = c("CTpan"),
  ED  = c("IT", "LE")
)

cat("  Zone-to-subcompartment mapping:\n")
for (sc in names(zone_to_subcomp)) {
  cat(sprintf("    %s <- %s\n", sc, paste(zone_to_subcomp[[sc]], collapse = " + ")))
}

# Transpose ssGSEA: samples x pathways
ssgsea_t <- as.data.frame(t(ssgsea_scores))
ssgsea_t$sample_id <- rownames(ssgsea_t)

# Add zone and patient_id from metadata
ssgsea_t$zone <- meta_matched$zone[match(ssgsea_t$sample_id,
                                          rownames(meta_matched))]
ssgsea_t$patient_id <- meta_matched$patient_id[match(ssgsea_t$sample_id,
                                                      rownames(meta_matched))]

# Handle case where match failed (try index-based)
if (all(is.na(ssgsea_t$zone))) {
  cat("  Direct sample matching failed. Using positional matching.\n")
  ssgsea_t$zone       <- meta_matched$zone
  ssgsea_t$patient_id <- meta_matched$patient_id
}

# Remove samples without zone/patient annotation
ssgsea_t <- ssgsea_t[!is.na(ssgsea_t$zone) & !is.na(ssgsea_t$patient_id), ]
cat(sprintf("  Samples with zone + patient annotation: %d\n", nrow(ssgsea_t)))

# Pathway column names
pathway_names <- rownames(ssgsea_scores)

# Aggregate: mean score per patient per subcompartment
agg_list <- list()

for (sc in names(zone_to_subcomp)) {
  zones <- zone_to_subcomp[[sc]]
  sc_data <- ssgsea_t[ssgsea_t$zone %in% zones, ]

  if (nrow(sc_data) == 0) {
    cat(sprintf("  WARNING: No samples for subcompartment %s (zones: %s)\n",
                sc, paste(zones, collapse = ", ")))
    next
  }

  # Mean across samples within each patient for this subcompartment
  sc_agg <- sc_data %>%
    group_by(patient_id) %>%
    summarise(across(all_of(pathway_names), ~ mean(.x, na.rm = TRUE)),
              n_samples = n(),
              .groups = "drop") %>%
    mutate(subcompartment = sc)

  cat(sprintf("  %s: %d patients (from %d samples, median %.0f samples/patient)\n",
              sc, nrow(sc_agg), nrow(sc_data), median(sc_agg$n_samples)))

  agg_list[[sc]] <- sc_agg
}

# Combine all subcompartments
ssgsea_aggregated <- bind_rows(agg_list)
cat(sprintf("  Aggregated ssGSEA: %d rows (patient x subcompartment)\n",
            nrow(ssgsea_aggregated)))

# Report patients with missing subcompartments
patients_per_sc <- ssgsea_aggregated %>%
  group_by(patient_id) %>%
  summarise(n_sc = n(), .groups = "drop")

cat(sprintf("  Patients with all 3 subcompartments: %d/%d\n",
            sum(patients_per_sc$n_sc == 3),
            length(unique(ssgsea_aggregated$patient_id))))
cat(sprintf("  Patients with 2 subcompartments: %d\n",
            sum(patients_per_sc$n_sc == 2)))
cat(sprintf("  Patients with 1 subcompartment: %d\n",
            sum(patients_per_sc$n_sc == 1)))

# ============================================================
# SECTION 2: PREPARE RADIOMIC FEATURES
# ============================================================

cat("\n--- Section 2: Prepare Radiomic Features ---\n")

# Load structured radiomics data (per-subcompartment, from script 00)
radio_per_sc <- readRDS(file.path(PROC_DIR, "radiomics_per_subcompartment.rds"))
cat(sprintf("  Loaded radiomics: %d subcompartments\n", length(radio_per_sc)))

# Deduplicate columns: bin-invariant features (Intensity, Morphologic, etc.)
# are repeated across bin configurations in Combined.xlsx — keep first occurrence
for (sc in names(radio_per_sc)) {
  cn <- colnames(radio_per_sc[[sc]])
  dupe_mask <- duplicated(cn)
  if (any(dupe_mask)) {
    cat(sprintf("  %s: removing %d duplicate columns (%d -> %d unique)\n",
                sc, sum(dupe_mask), length(cn), sum(!dupe_mask)))
    radio_per_sc[[sc]] <- radio_per_sc[[sc]][, !dupe_mask, drop = FALSE]
  }
}

# Normalize feature names: strip subcompartment identifier so columns align
# e.g., T1_T1_ET_Intensity_Energy -> T1_T1_Intensity_Energy (same across ET/NET/ED)
for (sc in names(radio_per_sc)) {
  cn <- colnames(radio_per_sc[[sc]])
  cn <- sub(paste0("_", sc, "_"), "_", cn)
  colnames(radio_per_sc[[sc]]) <- cn
}

# Verify alignment across subcompartments
ref_cols <- sort(setdiff(colnames(radio_per_sc[[1]]), "patient_id"))
for (sc in names(radio_per_sc)[-1]) {
  sc_cols <- sort(setdiff(colnames(radio_per_sc[[sc]]), "patient_id"))
  if (!identical(ref_cols, sc_cols)) {
    cat(sprintf("  WARNING: %s has %d/%d matching feature names with %s\n",
                sc, length(intersect(ref_cols, sc_cols)), length(ref_cols),
                names(radio_per_sc)[1]))
  }
}
cat(sprintf("  Feature names normalized across subcompartments: %d features\n",
            length(ref_cols)))

# Build long-format: stack subcompartments, subset to matched patients
feature_base_names <- setdiff(colnames(radio_per_sc[[1]]), "patient_id")

radio_long_list <- list()
for (sc in names(radio_per_sc)) {
  sc_df <- radio_per_sc[[sc]]
  sc_df$subcompartment <- sc
  # Subset to matched patients
  sc_df <- sc_df[sc_df$patient_id %in% matched_patients, ]
  radio_long_list[[sc]] <- sc_df
  cat(sprintf("  %s: %d matched patients, %d features\n",
              sc, nrow(sc_df), length(feature_base_names)))
}

radio_long <- bind_rows(radio_long_list)
cat(sprintf("  Radiomics long format (matched): %d rows x %d feature cols\n",
            nrow(radio_long), length(feature_base_names)))
cat(sprintf("  Patients: %d, Subcompartments: %s\n",
            length(unique(radio_long$patient_id)),
            paste(sort(unique(radio_long$subcompartment)), collapse = ", ")))

# ============================================================
# SECTION 3: MERGE AND FORMAT FOR SCRIPT 03
# ============================================================

cat("\n--- Section 3: Merge and Format ---\n")

# Rename pathway columns to pw_1 .. pw_24
pw_original <- pathway_names <- setdiff(
  colnames(ssgsea_aggregated),
  c("patient_id", "subcompartment", "n_samples")
)

pw_new <- paste0("pw_", seq_along(pw_original))
pw_lookup <- data.frame(
  pw_col = pw_new,
  pathway_name = pw_original,
  stringsAsFactors = FALSE
)
cat("  Pathway column mapping:\n")
for (i in seq_along(pw_original)) {
  cat(sprintf("    %s -> %s\n", pw_new[i], pw_original[i]))
}

# Rename in ssgsea data
ssgsea_for_merge <- ssgsea_aggregated %>%
  select(patient_id, subcompartment, all_of(pw_original))
colnames(ssgsea_for_merge) <- c("patient_id", "subcompartment", pw_new)

# Rename radiomic features to rad_1 .. rad_N
rad_original <- feature_base_names
rad_new <- paste0("rad_", seq_along(rad_original))
rad_lookup <- data.frame(
  rad_col = rad_new,
  feature_name = rad_original,
  stringsAsFactors = FALSE
)

# Rename in radiomics data
radio_for_merge <- radio_long %>%
  select(patient_id, subcompartment, all_of(rad_original))
colnames(radio_for_merge) <- c("patient_id", "subcompartment", rad_new)

cat(sprintf("  Pathway columns: %d (pw_1 to pw_%d)\n",
            length(pw_new), length(pw_new)))
cat(sprintf("  Radiomic columns: %d (rad_1 to rad_%d)\n",
            length(rad_new), length(rad_new)))

# Inner join on patient_id + subcompartment
analysis_data <- inner_join(
  ssgsea_for_merge,
  radio_for_merge,
  by = c("patient_id", "subcompartment")
)

cat(sprintf("  Merged dataset: %d rows x %d columns\n",
            nrow(analysis_data), ncol(analysis_data)))
cat(sprintf("  Patients in merged data: %d\n",
            length(unique(analysis_data$patient_id))))

# Set subcompartment as factor with ET as reference
analysis_data$subcompartment <- factor(
  analysis_data$subcompartment,
  levels = c("ET", "NET", "ED")
)

# Report rows per subcompartment
cat("  Rows per subcompartment:\n")
for (sc in levels(analysis_data$subcompartment)) {
  n <- sum(analysis_data$subcompartment == sc)
  cat(sprintf("    %s: %d\n", sc, n))
}

# ============================================================
# SECTION 4: FINAL VALIDATION
# ============================================================

cat("\n--- Section 4: Final Validation ---\n")

# Match script 03's assertions (lines 95-101)
pathway_cols_check  <- grep("^pw_", colnames(analysis_data), value = TRUE)
radiomic_cols_check <- grep("^rad_", colnames(analysis_data), value = TRUE)

cat(sprintf("  patient_id column: %s\n",
            ifelse("patient_id" %in% colnames(analysis_data), "OK", "MISSING")))
cat(sprintf("  subcompartment column: %s\n",
            ifelse("subcompartment" %in% colnames(analysis_data), "OK", "MISSING")))
cat(sprintf("  Pathway columns (pw_*): %d (expected 24)\n",
            length(pathway_cols_check)))
cat(sprintf("  Radiomic columns (rad_*): %d\n",
            length(radiomic_cols_check)))
cat(sprintf("  Unique patients: %d (expected >= 25)\n",
            length(unique(analysis_data$patient_id))))

# Run script 03's stopifnot checks
stopifnot(
  "patient_id column missing" =
    "patient_id" %in% colnames(analysis_data),
  "subcompartment column missing" =
    "subcompartment" %in% colnames(analysis_data),
  "Expected 24 pathway columns" =
    length(pathway_cols_check) == 24,
  "No radiomic columns found" =
    length(radiomic_cols_check) > 0,
  "Expected >= 25 patients" =
    length(unique(analysis_data$patient_id)) >= 25
)
cat("  All script 03 assertions PASSED.\n")

# Check for NAs in critical columns
pw_na <- sum(is.na(analysis_data[, pathway_cols_check]))
rad_na <- sum(is.na(analysis_data[, radiomic_cols_check]))
cat(sprintf("  NA values in pathway columns: %d\n", pw_na))
cat(sprintf("  NA values in radiomic columns: %d\n", rad_na))

if (rad_na > 0) {
  pct_na <- 100 * rad_na / (nrow(analysis_data) * length(radiomic_cols_check))
  cat(sprintf("  Radiomic NA percentage: %.2f%%\n", pct_na))
  if (pct_na > 10) {
    cat("  WARNING: >10%% NAs in radiomic features. Consider imputation.\n")
  }
}

# Replace Inf with NA
n_inf <- sum(is.infinite(as.matrix(analysis_data[, radiomic_cols_check])))
if (n_inf > 0) {
  cat(sprintf("  Replacing %d infinite values with NA.\n", n_inf))
  for (col in radiomic_cols_check) {
    analysis_data[[col]][is.infinite(analysis_data[[col]])] <- NA
  }
}

# ============================================================
# SAVE OUTPUTS
# ============================================================

cat("\n--- Saving Outputs ---\n")

# Save analysis-ready dataset
output_file <- file.path(DATA_DIR, "analysis_data.rds")
saveRDS(analysis_data, output_file)
cat(sprintf("  analysis_data.rds saved: %s\n", output_file))
cat(sprintf("  Dimensions: %d rows x %d columns\n",
            nrow(analysis_data), ncol(analysis_data)))

# Save lookup tables
saveRDS(pw_lookup, file.path(PROC_DIR, "pathway_lookup.rds"))
saveRDS(rad_lookup, file.path(PROC_DIR, "feature_lookup.rds"))
cat("  Lookup tables saved (pathway_lookup.rds, feature_lookup.rds)\n")

# Print summary for verification
cat("\n  === FINAL SUMMARY ===\n")
cat(sprintf("  Patients: %d\n", length(unique(analysis_data$patient_id))))
cat(sprintf("  Rows (patient x subcompartment): %d\n", nrow(analysis_data)))
cat(sprintf("  Pathway scores: %d (pw_1 to pw_%d)\n",
            length(pathway_cols_check), length(pathway_cols_check)))
cat(sprintf("  Radiomic features: %d (rad_1 to rad_%d)\n",
            length(radiomic_cols_check), length(radiomic_cols_check)))
cat(sprintf("  Subcompartment levels: %s\n",
            paste(levels(analysis_data$subcompartment), collapse = ", ")))
cat(sprintf("  Ready for: Rscript R-RADIOMAP-IvyGAP/03_statistical_analysis.R\n"))

cat("\n============================================================\n")
cat("Script 02 complete. analysis_data.rds ready for script 03.\n")
cat("============================================================\n")
