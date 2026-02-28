# ============================================================
# RADIOMAP-IvyGAP: ssGSEA Pathway Enrichment
# Script 01: Compute 24 ssGSEA scores for matched IvyGAP samples
# ============================================================
# Author: Daniele Piccolo, MD
# Date: 2026-02-21
# Dependencies: GSVA, msigdbr, SummarizedExperiment, dplyr, ggplot2,
#               tidyr, ivygapSE
# Input:
#   data/raw/rnaseq/ivygap_expression_matrix.rds
#   data/raw/rnaseq/ivygap_sample_metadata.rds
#   data/processed/matched_patients.csv
#   data/processed/audit_variables.rds
# Output:
#   data/processed/ssgsea_scores.rds (24 x N_samples matrix)
#   data/processed/gene_sets_used.rds (named list for Jaccard overlap)
#   figures/ssgsea_qc_zone_modules.png
#   figures/ssgsea_qc_zone_modules.pdf
# ============================================================

cat("============================================================\n")
cat("RADIOMAP-IvyGAP: Script 01 — ssGSEA Computation\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("============================================================\n\n")

set.seed(42)

# --- Setup ---
PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RAW_DIR     <- file.path(DATA_DIR, "raw")
PROC_DIR    <- file.path(DATA_DIR, "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Install/load packages
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

install_if_missing("GSVA", bioc = TRUE)
install_if_missing("msigdbr")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("ggplot2")

cat("  All packages loaded.\n\n")

# ============================================================
# SECTION 1: PREPARE EXPRESSION MATRIX
# ============================================================

cat("--- Section 1: Prepare Expression Matrix ---\n")

# Load script 00 outputs
expr_matrix <- readRDS(file.path(RAW_DIR, "rnaseq", "ivygap_expression_matrix.rds"))
sample_meta <- readRDS(file.path(RAW_DIR, "rnaseq", "ivygap_sample_metadata.rds"))
audit_vars  <- readRDS(file.path(PROC_DIR, "audit_variables.rds"))

matched_patients <- audit_vars$matched_patients
zone_col         <- audit_vars$zone_col
primary_zones    <- audit_vars$primary_zones

cat(sprintf("  Full expression matrix: %d genes x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# Normalize zone labels (same as script 00)
sample_meta$zone <- sapply(strsplit(as.character(sample_meta[[zone_col]]), "-"), `[`, 1)

# Subset to matched patients
sample_meta$patient_id <- sub("^(W\\d+).*", "\\1",
                              as.character(sample_meta[[audit_vars$tumor_col]]))
matched_idx <- sample_meta$patient_id %in% matched_patients
sample_meta_matched <- sample_meta[matched_idx, ]

# Match expression matrix columns to metadata rows
# Column names in expr_matrix should correspond to sample IDs in metadata
common_samples <- intersect(colnames(expr_matrix), rownames(sample_meta_matched))
if (length(common_samples) == 0) {
  # Try matching on a different identifier
  cat("  No direct column-rowname match. Trying index-based matching...\n")
  # Assume same order if same length
  if (ncol(expr_matrix) == nrow(sample_meta)) {
    expr_matched <- expr_matrix[, matched_idx]
    meta_matched <- sample_meta_matched
  } else {
    stop("Cannot match expression columns to metadata rows. ",
         "Expression has ", ncol(expr_matrix), " columns; ",
         "metadata has ", nrow(sample_meta), " rows.")
  }
} else {
  expr_matched <- expr_matrix[, common_samples]
  meta_matched <- sample_meta_matched[common_samples, ]
}

cat(sprintf("  Matched expression matrix: %d genes x %d samples\n",
            nrow(expr_matched), ncol(expr_matched)))
cat(sprintf("  Patients represented: %d\n",
            length(unique(meta_matched$patient_id))))
cat(sprintf("  Zones represented: %s\n",
            paste(sort(unique(meta_matched$zone)), collapse = ", ")))

# Check gene identifiers (symbols vs Entrez IDs vs Ensembl)
sample_genes <- head(rownames(expr_matched), 10)
cat(sprintf("  Gene ID format (first 10): %s\n",
            paste(sample_genes, collapse = ", ")))

# Determine if gene IDs are symbols, Entrez, or Ensembl
is_symbol  <- any(grepl("^[A-Z][A-Z0-9]+$", sample_genes))
is_entrez  <- all(grepl("^\\d+$", sample_genes))
is_ensembl <- any(grepl("^ENSG", sample_genes))

if (is_entrez) {
  cat("  Gene IDs appear to be Entrez IDs. Converting to symbols...\n")
  # Try to load org.Hs.eg.db for conversion
  install_if_missing("org.Hs.eg.db", bioc = TRUE)
  install_if_missing("AnnotationDbi", bioc = TRUE)
  entrez_to_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expr_matched),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  # Remove NAs and duplicates
  valid <- !is.na(entrez_to_symbol) & !duplicated(entrez_to_symbol)
  expr_matched <- expr_matched[valid, ]
  rownames(expr_matched) <- entrez_to_symbol[valid]
  cat(sprintf("  After conversion: %d genes with symbols\n", nrow(expr_matched)))
} else if (is_ensembl) {
  cat("  Gene IDs appear to be Ensembl IDs. Converting to symbols...\n")
  install_if_missing("org.Hs.eg.db", bioc = TRUE)
  install_if_missing("AnnotationDbi", bioc = TRUE)
  ensembl_to_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = sub("\\..*", "", rownames(expr_matched)),  # remove version
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  valid <- !is.na(ensembl_to_symbol) & !duplicated(ensembl_to_symbol)
  expr_matched <- expr_matched[valid, ]
  rownames(expr_matched) <- ensembl_to_symbol[valid]
  cat(sprintf("  After conversion: %d genes with symbols\n", nrow(expr_matched)))
} else {
  cat("  Gene IDs appear to be gene symbols. No conversion needed.\n")
}

# Filter low-expression genes: FPKM > 1 in >= 10% of samples
min_samples <- ceiling(0.10 * ncol(expr_matched))
gene_pass <- rowSums(expr_matched > 1) >= min_samples
expr_filtered <- expr_matched[gene_pass, ]
cat(sprintf("  Gene filter (FPKM>1 in >=10%% samples): %d -> %d genes\n",
            nrow(expr_matched), nrow(expr_filtered)))

# Log2 transform for ssGSEA with kcdf="Gaussian"
expr_log2 <- log2(expr_filtered + 1)
cat(sprintf("  Transformed: log2(FPKM + 1), range [%.2f, %.2f]\n",
            min(expr_log2), max(expr_log2)))
cat("\n")

# ============================================================
# SECTION 2: ASSEMBLE 24 GENE SETS
# ============================================================

cat("--- Section 2: Assemble 24 Gene Sets ---\n")

gene_sets <- list()

# --- 2a. 15 Hallmark gene sets from MSigDB ---
cat("  Loading 15 Hallmark gene sets from MSigDB...\n")
hallmark_db <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")

hallmark_names <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_COMPLEMENT"
)

# Short display names
hallmark_short <- c(
  "Hypoxia", "Angiogenesis", "EMT", "Inflammatory_Response",
  "TNFA_NFkB", "IL6_JAK_STAT3", "IFN_Gamma", "P53_Pathway",
  "MYC_V1", "E2F_Targets", "G2M_Checkpoint", "mTORC1",
  "Glycolysis", "OxPhos", "Complement"
)

for (i in seq_along(hallmark_names)) {
  genes <- hallmark_db$gene_symbol[hallmark_db$gs_name == hallmark_names[i]]
  genes <- unique(genes)
  gene_sets[[hallmark_short[i]]] <- genes
  cat(sprintf("    %s: %d genes\n", hallmark_short[i], length(genes)))
}

# --- 2b. 4 Neftel 2019 cellular state signatures ---
# Hard-coded from Neftel et al. 2019, Cell 178(4):835-849, Table S5
# MES1+MES2 collapsed into MES; NPC1+NPC2 collapsed into NPC
cat("  Assembling 4 Neftel cellular state signatures...\n")

# MES-like (MES1 + MES2 meta-modules)
neftel_MES <- c(
  # MES1
  "CHI3L1", "VIM", "ANXA2", "CD44", "LGALS3", "S100A11", "TGFBI",
  "THBS1", "SERPINE1", "ADM", "PDPN", "TAGLN2", "IGFBP3", "MGP",
  "CLIC1", "MYL9", "ACTN1", "CALD1", "IGFBP5", "DCBLD2", "FLNC",
  "S100A10", "ANXA1", "COL1A2", "COL5A1", "COL3A1",
  # MES2
  "POSTN", "TNC", "FN1", "COL1A1", "COL4A1", "COL4A2", "COL6A2",
  "COL6A1", "SPARC", "LOXL2", "PLOD2", "ITGA5", "MMP2", "LOX",
  "PMEPA1", "NNMT", "BCL3", "NAMPT", "IL6", "TGFB1", "TNFAIP6",
  "DAB2", "TIMP1", "EMP1", "ANXA5", "S100A6", "S100A4", "CTGF"
)

# AC-like (astrocyte-like meta-module)
neftel_AC <- c(
  "GFAP", "SLC1A3", "AQP4", "ALDOC", "CLU", "MT1X", "MT2A",
  "SOX9", "NTRK2", "SLC1A2", "GJA1", "FGFR3", "ID4", "APOE",
  "S100B", "GLUL", "SLC4A4", "ATP1A2", "AGT", "BCAN", "HEY1",
  "CST3", "ATP1B2", "PLPP3", "MERTK", "DTNA", "BBOX1", "SOX2",
  "PRCP", "PON2", "TJP1", "GRAMD3", "DBI", "SLC7A11", "NCAN",
  "SPARCL1", "GPM6A", "F3", "NDRG2", "ALDH1L1", "SOX21", "FABP7",
  "HEPN1", "CPE", "CABLES1", "PAPLN", "ETNPPL", "PRODH", "ACSL6",
  "GRM3"
)

# OPC-like (oligodendrocyte precursor meta-module)
neftel_OPC <- c(
  "PDGFRA", "OLIG1", "OLIG2", "SOX10", "NKX2-2", "CSPG4", "GPR17",
  "PLP1", "MBP", "CNP", "TNR", "LHFPL3", "PCDH15", "TMEM2",
  "NEU4", "LUZP2", "EPN2", "OPCML", "NXPH1", "RAB33A", "OMG",
  "BCAS1", "SEMA5A", "PLLP", "EDIL3", "SCRG1", "DLL3", "ASCL1",
  "SOX8", "HIP1R", "FXYD6", "BCHE", "APOD", "FEZ1", "CDH20",
  "FIBIN", "TTYH1", "GNG8", "SNX22", "FERMT1", "CA10", "GALNT13",
  "LIMA1", "LPPR1", "SIRT2", "NFASC", "FA2H", "CLDN11", "UGT8",
  "MAG"
)

# NPC-like (NPC1 + NPC2 neural progenitor meta-modules)
neftel_NPC <- c(
  # NPC1
  "DCX", "DLX1", "DLX2", "DLX5", "DLX6", "SOX4", "SOX11",
  "STMN2", "TUBB3", "MAP2", "NCAM1", "CD24", "NRXN3", "GAD1",
  "GAD2", "SYT1", "SNAP25", "NSG1", "STMN1", "MLLT11", "DPYSL3",
  "ST18", "NHLH1", "INSM1", "HES6",
  # NPC2
  "PROX1", "SP8", "ARX", "SRRM4", "CELF4", "PPP2R2B", "NOVA1",
  "ELAVL3", "ELAVL4", "RBFOX1", "MEIS2", "PBX3", "ZNF711",
  "CNTN2", "NRCAM", "SCG2", "SVOP", "SYP", "ATP1A3", "CAMK2B",
  "NEUROD1", "NEUROD2", "NEUROD6", "SOX5", "LHX2", "EMX1",
  "TBR1", "SATB2", "BCL11B"
)

gene_sets[["Neftel_MES"]] <- unique(neftel_MES)
gene_sets[["Neftel_AC"]]  <- unique(neftel_AC)
gene_sets[["Neftel_OPC"]] <- unique(neftel_OPC)
gene_sets[["Neftel_NPC"]] <- unique(neftel_NPC)

cat(sprintf("    Neftel_MES: %d genes\n", length(gene_sets[["Neftel_MES"]])))
cat(sprintf("    Neftel_AC: %d genes\n", length(gene_sets[["Neftel_AC"]])))
cat(sprintf("    Neftel_OPC: %d genes\n", length(gene_sets[["Neftel_OPC"]])))
cat(sprintf("    Neftel_NPC: %d genes\n", length(gene_sets[["Neftel_NPC"]])))

# --- 2c. 5 IvyGAP zone-specific modules (data-derived) ---
# Top 200 DE genes per zone vs all others (Wilcoxon, FDR < 0.05, |log2FC| > 1)
cat("  Deriving 5 IvyGAP zone-specific gene modules...\n")

# Use ALL samples (not just matched) for DE to maximize power
sample_meta_all <- sample_meta
sample_meta_all$zone <- sapply(strsplit(as.character(sample_meta_all[[zone_col]]), "-"), `[`, 1)

# Match samples to expression columns
if (all(rownames(sample_meta_all) %in% colnames(expr_matrix))) {
  expr_for_de <- expr_matrix[, rownames(sample_meta_all)]
  meta_for_de <- sample_meta_all
} else {
  # Use positional matching
  expr_for_de <- expr_matrix
  meta_for_de <- sample_meta_all
}

# Apply same gene filter
gene_pass_de <- rowSums(expr_for_de > 1) >= ceiling(0.10 * ncol(expr_for_de))
expr_for_de <- expr_for_de[gene_pass_de, ]
# Log2 transform for fold change
expr_log2_de <- log2(expr_for_de + 1)

zone_modules <- list()
de_zones <- c("CT", "CTmvp", "CTpan", "IT", "LE")

for (z in de_zones) {
  cat(sprintf("    Computing DE for zone %s...\n", z))

  is_zone <- meta_for_de$zone == z
  n_in  <- sum(is_zone)
  n_out <- sum(!is_zone)

  if (n_in < 5) {
    cat(sprintf("      WARNING: Only %d samples in zone %s. Skipping.\n", n_in, z))
    next
  }

  # Wilcoxon rank-sum test per gene
  pvals <- numeric(nrow(expr_log2_de))
  log2fc <- numeric(nrow(expr_log2_de))

  for (g in seq_len(nrow(expr_log2_de))) {
    in_vals  <- expr_log2_de[g, is_zone]
    out_vals <- expr_log2_de[g, !is_zone]

    log2fc[g] <- mean(in_vals) - mean(out_vals)

    wt <- suppressWarnings(
      wilcox.test(in_vals, out_vals, exact = FALSE)
    )
    pvals[g] <- wt$p.value
  }

  # FDR correction
  fdr <- p.adjust(pvals, method = "BH")

  # Filter: FDR < 0.05 AND |log2FC| > 1 AND upregulated (positive FC)
  sig_up <- which(fdr < 0.05 & log2fc > 1)

  if (length(sig_up) == 0) {
    # Relax to FDR < 0.10
    sig_up <- which(fdr < 0.10 & log2fc > 0.5)
    cat(sprintf("      Relaxed criteria (FDR<0.10, |FC|>0.5): %d genes\n",
                length(sig_up)))
  }

  if (length(sig_up) > 0) {
    # Sort by fold change descending, take top 200
    sig_up_ordered <- sig_up[order(log2fc[sig_up], decreasing = TRUE)]
    top_genes <- rownames(expr_log2_de)[head(sig_up_ordered, 200)]
    zone_modules[[z]] <- top_genes
    cat(sprintf("      %s: %d significant genes, top %d selected\n",
                z, length(sig_up), length(top_genes)))
  } else {
    cat(sprintf("      WARNING: No significant DE genes for zone %s.\n", z))
    zone_modules[[z]] <- character(0)
  }
}

# Add zone modules to gene sets with descriptive names
zone_module_names <- c(
  CT     = "IvyGAP_CT_module",
  CTmvp  = "IvyGAP_CTmvp_module",
  CTpan  = "IvyGAP_CTpan_module",
  IT     = "IvyGAP_IT_module",
  LE     = "IvyGAP_LE_module"
)

for (z in names(zone_modules)) {
  gs_name <- zone_module_names[z]
  gene_sets[[gs_name]] <- zone_modules[[z]]
  cat(sprintf("    %s: %d genes\n", gs_name, length(zone_modules[[z]])))
}

# Report gene set summary
cat(sprintf("\n  Total gene sets assembled: %d\n", length(gene_sets)))
cat(sprintf("    Hallmark: %d\n", sum(names(gene_sets) %in% hallmark_short)))
cat(sprintf("    Neftel: %d\n",
            sum(grepl("^Neftel_", names(gene_sets)))))
cat(sprintf("    IvyGAP zone: %d\n",
            sum(grepl("^IvyGAP_", names(gene_sets)))))

# Check gene overlap with expression matrix
for (gs_name in names(gene_sets)) {
  n_total <- length(gene_sets[[gs_name]])
  n_found <- sum(gene_sets[[gs_name]] %in% rownames(expr_log2))
  if (n_found < 10) {
    cat(sprintf("  WARNING: %s — only %d/%d genes found in expression data\n",
                gs_name, n_found, n_total))
  }
}

stopifnot(
  "Expected 24 gene sets" = length(gene_sets) == 24
)
cat("  VALIDATED: 24 gene sets assembled.\n\n")

# ============================================================
# SECTION 3: RUN ssGSEA
# ============================================================

cat("--- Section 3: Run ssGSEA ---\n")

# Use matched samples only for ssGSEA (these go into analysis)
cat(sprintf("  Running ssGSEA on %d genes x %d samples...\n",
            nrow(expr_log2), ncol(expr_log2)))
cat("  Parameters: method=ssgsea, kcdf=Gaussian, min.sz=10, ssgsea.norm=TRUE\n")
cat("  This may take a few minutes...\n")

set.seed(42)

# GSVA::gsva() interface (handles both old and new API)
ssgsea_result <- tryCatch({
  # New GSVA API (>= 1.50.0): uses gsvaParam object
  param <- GSVA::ssgseaParam(
    exprData = expr_log2,
    geneSets = gene_sets,
    minSize  = 10,
    normalize = TRUE
  )
  GSVA::gsva(param)
}, error = function(e) {
  cat(sprintf("  New API failed (%s), trying legacy API...\n", e$message))
  # Legacy GSVA API (< 1.50.0)
  GSVA::gsva(
    expr       = expr_log2,
    gset.idx.list = gene_sets,
    method     = "ssgsea",
    kcdf       = "Gaussian",
    min.sz     = 10,
    ssgsea.norm = TRUE,
    verbose    = TRUE
  )
})

cat(sprintf("  ssGSEA complete: %d pathways x %d samples\n",
            nrow(ssgsea_result), ncol(ssgsea_result)))
cat(sprintf("  Pathways computed: %s\n",
            paste(rownames(ssgsea_result), collapse = ", ")))

# Check for NAs
n_na_ssgsea <- sum(is.na(ssgsea_result))
if (n_na_ssgsea > 0) {
  cat(sprintf("  WARNING: %d NA values in ssGSEA scores.\n", n_na_ssgsea))
}
cat("\n")

# ============================================================
# SECTION 4: QC — ZONE MODULE VALIDATION
# ============================================================

cat("--- Section 4: QC Plots ---\n")

# Zone modules should show highest enrichment in their own zones
# This serves as internal validation of the ssGSEA computation

# Build data frame for plotting
zone_module_rows <- grep("^IvyGAP_", rownames(ssgsea_result), value = TRUE)

if (length(zone_module_rows) > 0) {
  ssgsea_df <- as.data.frame(t(ssgsea_result[zone_module_rows, , drop = FALSE]))
  ssgsea_df$sample_id <- rownames(ssgsea_df)

  # Merge with zone labels
  ssgsea_df$zone <- meta_matched$zone[match(ssgsea_df$sample_id,
                                             rownames(meta_matched))]

  # Handle case where match failed (try column names)
  if (all(is.na(ssgsea_df$zone))) {
    # Try index-based
    ssgsea_df$zone <- meta_matched$zone
  }

  # Remove samples without zone annotation
  ssgsea_df <- ssgsea_df[!is.na(ssgsea_df$zone), ]

  # Pivot to long format
  ssgsea_long <- ssgsea_df %>%
    pivot_longer(
      cols      = all_of(zone_module_rows),
      names_to  = "module",
      values_to = "score"
    ) %>%
    mutate(
      module_zone = sub("IvyGAP_(.+)_module", "\\1", module),
      is_own_zone = (zone == module_zone)
    )

  # Boxplot: scores by zone, faceted by module
  p_qc <- ggplot(ssgsea_long, aes(x = zone, y = score, fill = is_own_zone)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    facet_wrap(~ module, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c("FALSE" = "gray80", "TRUE" = "steelblue"),
      labels = c("Other zones", "Own zone"),
      name   = ""
    ) +
    labs(
      title    = "IvyGAP Zone Module Validation",
      subtitle = "Zone-specific modules should show highest scores in their own zones",
      x = "IvyGAP Zone",
      y = "ssGSEA Enrichment Score"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
      strip.text     = element_text(size = 8),
      legend.position = "bottom"
    )

  ggsave(file.path(FIG_DIR, "ssgsea_qc_zone_modules.png"), p_qc,
         width = 10, height = 8, dpi = 300)
  ggsave(file.path(FIG_DIR, "ssgsea_qc_zone_modules.pdf"), p_qc,
         width = 10, height = 8)
  cat("  QC plots saved to figures/\n")

  # Report: does each module discriminate its own zone?
  cat("  Zone module discrimination (own zone vs others):\n")
  for (mod in zone_module_rows) {
    mod_zone <- sub("IvyGAP_(.+)_module", "\\1", mod)
    own  <- ssgsea_long$score[ssgsea_long$module == mod & ssgsea_long$is_own_zone]
    other <- ssgsea_long$score[ssgsea_long$module == mod & !ssgsea_long$is_own_zone]
    if (length(own) > 2 && length(other) > 2) {
      wt <- wilcox.test(own, other, alternative = "greater")
      effect <- median(own) - median(other)
      cat(sprintf("    %s: median diff = %.3f, Wilcox p = %.2e %s\n",
                  mod, effect, wt$p.value,
                  ifelse(wt$p.value < 0.05, "(PASS)", "(FAIL)")))
    }
  }
} else {
  cat("  WARNING: No IvyGAP zone modules found in ssGSEA results.\n")
}

# ============================================================
# SAVE OUTPUTS
# ============================================================

cat("\n--- Saving Outputs ---\n")

# Save ssGSEA scores matrix
saveRDS(ssgsea_result, file.path(PROC_DIR, "ssgsea_scores.rds"))
cat(sprintf("  ssGSEA scores saved: %s\n",
            file.path(PROC_DIR, "ssgsea_scores.rds")))

# Save gene sets (for Jaccard overlap analysis in script 03)
saveRDS(gene_sets, file.path(PROC_DIR, "gene_sets_used.rds"))
cat(sprintf("  Gene sets saved: %s\n",
            file.path(PROC_DIR, "gene_sets_used.rds")))

# Save metadata for matched samples (needed by script 02)
saveRDS(meta_matched, file.path(PROC_DIR, "sample_metadata_matched.rds"))
cat(sprintf("  Matched sample metadata saved.\n"))

cat("\n============================================================\n")
cat(sprintf("Script 01 complete. %d pathways x %d samples.\n",
            nrow(ssgsea_result), ncol(ssgsea_result)))
cat("============================================================\n")
