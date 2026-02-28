# ============================================================
# RADIOMAP-IvyGAP: Nested Cross-Validation Analysis
# Script 03b: Bias-Free Predictive Analysis + Clinical Covariates
# ============================================================
# Author: Daniele Piccolo, MD
# Date: 2026-02-22
# Purpose: Address data leakage concern — the original pipeline
#   performed univariate feature screening on the full dataset
#   before LOPO-CV, inflating apparent associations. Here, feature
#   selection is performed INSIDE each CV fold (nested CV), making
#   this the PRIMARY predictive analysis.
#
# Analyses:
#   1. Nested LOPO-CV with internal feature selection (all 24 pathways)
#   2. Feature stability analysis (selection frequency across folds)
#   3. Bootstrap CIs for nested CV performance
#   4. Nested permutation test (unbiased p-values for R2_cv > 0 pathways)
#   5. Clinical covariate adjustment (S9: age, MGMT, molecular subtype)
#   6. Feature lookup table export for supplementary materials
#   7. Comparison table: nested CV vs. legacy pre-screened results
#
# Input:  results/03_analysis_session.RData
#         data/processed/feature_lookup.rds
#         data/processed/pathway_lookup.rds
#         data/raw/rnaseq/ivygap_tumor_details.rds
# Output: results/nested_cv_results.rds
#         results/nested_cv_summary.csv
#         results/nested_permutation_results.rds
#         results/clinical_covariates_s9.rds
#         results/feature_lookup.csv
#         results/pw4_feature_stability.csv
#         results/nested_vs_legacy_comparison.csv
#         results/03b_nested_cv_session.RData
# ============================================================

# --- Setup ---
library(lme4)
library(lmerTest)
library(performance)
library(glmnet)
library(dplyr)
library(tidyr)
library(broom.mixed)

set.seed(42)

PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

cat("============================================================\n")
cat("RADIOMAP-IvyGAP: Script 03b - Nested CV Analysis\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("============================================================\n\n")

# ============================================================
# STEP 0: LOAD SESSION DATA
# ============================================================

cat("========== STEP 0: LOAD SESSION DATA ==========\n")

load(file.path(RESULTS_DIR, "03_analysis_session.RData"))

cat(sprintf("Loaded session: %d patients, %d pathways, %d features (post-cor-filter)\n",
            length(unique(analysis_data$patient_id)),
            length(pathway_cols),
            length(features_pass_cor)))

# Load lookup tables
rad_lookup <- readRDS(file.path(DATA_DIR, "processed", "feature_lookup.rds"))
pw_lookup  <- readRDS(file.path(DATA_DIR, "processed", "pathway_lookup.rds"))

# ============================================================
# STEP 1: NESTED LOPO-CV WITH INTERNAL FEATURE SELECTION
# ============================================================

cat("\n========== STEP 1: NESTED LOPO-CV (ALL 24 PATHWAYS) ==========\n")

# Use patient-level data (already aggregated in script 03)
# z-scored features and pathways at patient level
all_z_features <- paste0("z_", features_pass_cor)
N_patients <- nrow(patient_level)

cat(sprintf("Patient-level data: N=%d, features=%d\n", N_patients, length(all_z_features)))

nested_cv_results <- list()

for (pw in pathway_cols) {
  z_pathway <- paste0("z_", pw)
  y <- patient_level[[z_pathway]]
  X_all <- as.matrix(patient_level[, all_z_features, drop = FALSE])
  N <- length(y)

  predictions <- numeric(N)
  fold_features <- vector("list", N)
  fold_alpha <- numeric(N)
  fold_lambda <- numeric(N)
  fold_n_sig <- integer(N)
  fold_coefs <- vector("list", N)

  for (i in 1:N) {
    X_train <- X_all[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X_all[i, , drop = FALSE]

    # Step A: Univariate Spearman filter on TRAINING set only
    spearman_p <- apply(X_train, 2, function(col) {
      tryCatch(cor.test(col, y_train, method = "spearman", exact = FALSE)$p.value,
               error = function(e) 1)
    })
    fdr_vals <- p.adjust(spearman_p, method = "BH")
    sig_feats <- names(fdr_vals)[fdr_vals < 0.10]

    # Take top 5 by raw p-value if > 5 pass
    if (length(sig_feats) > 5) {
      sig_feats <- names(sort(spearman_p[sig_feats]))[1:5]
    }

    fold_n_sig[i] <- length(sig_feats)
    fold_features[[i]] <- sig_feats

    if (length(sig_feats) == 0) {
      predictions[i] <- mean(y_train)
      fold_alpha[i] <- NA
      fold_lambda[i] <- NA
      fold_coefs[[i]] <- NULL
      next
    }

    # Step B: Elastic Net on selected features
    X_sel_train <- X_train[, sig_feats, drop = FALSE]
    X_sel_test  <- X_test[, sig_feats, drop = FALSE]

    best_alpha <- 0.5; best_lambda <- NULL; best_cvm <- Inf
    for (av in seq(0.1, 1.0, by = 0.1)) {
      cvf <- tryCatch(
        cv.glmnet(X_sel_train, y_train, alpha = av, nfolds = min(5, N - 1)),
        error = function(e) NULL)
      if (!is.null(cvf)) {
        idx_1se <- which(cvf$lambda == cvf$lambda.1se)
        if (length(idx_1se) > 0 && cvf$cvm[idx_1se] < best_cvm) {
          best_cvm <- cvf$cvm[idx_1se]; best_alpha <- av
          best_lambda <- cvf$lambda.1se
        }
      }
    }

    if (!is.null(best_lambda)) {
      fit <- glmnet(X_sel_train, y_train, alpha = best_alpha, lambda = best_lambda)
      predictions[i] <- predict(fit, newx = X_sel_test)[1]
      fold_alpha[i] <- best_alpha
      fold_lambda[i] <- best_lambda
      # Store non-zero coefficients
      cf <- as.matrix(coef(fit))
      fold_coefs[[i]] <- cf[cf[,1] != 0, , drop = FALSE]
    } else {
      predictions[i] <- mean(y_train)
      fold_alpha[i] <- NA
      fold_lambda[i] <- NA
      fold_coefs[[i]] <- NULL
    }
  }

  # Performance metrics
  SS_res <- sum((y - predictions)^2)
  SS_tot <- sum((y - mean(y))^2)
  R2_cv  <- 1 - SS_res / SS_tot
  MAE    <- mean(abs(y - predictions))
  rho    <- tryCatch(cor(y, predictions, method = "spearman"), error = function(e) NA)

  # Feature stability: count how often each feature was selected
  all_selected <- unlist(fold_features)
  if (length(all_selected) > 0) {
    feat_freq <- table(all_selected)
    feat_stability <- data.frame(
      feature   = names(feat_freq),
      n_folds   = as.integer(feat_freq),
      pct_folds = as.numeric(feat_freq) / N * 100,
      stringsAsFactors = FALSE
    ) %>% arrange(desc(n_folds))
  } else {
    feat_stability <- data.frame(
      feature = character(), n_folds = integer(),
      pct_folds = numeric(), stringsAsFactors = FALSE
    )
  }

  # Features selected in >50% of folds
  stable_features <- feat_stability$feature[feat_stability$pct_folds > 50]

  nested_cv_results[[pw]] <- list(
    R2_cv           = R2_cv,
    MAE             = MAE,
    spearman        = rho,
    predictions     = predictions,
    observed        = y,
    fold_features   = fold_features,
    fold_alpha      = fold_alpha,
    fold_lambda     = fold_lambda,
    fold_n_sig      = fold_n_sig,
    fold_coefs      = fold_coefs,
    feat_stability  = feat_stability,
    stable_features = stable_features
  )

  cat(sprintf("  %s: R2_cv=%.3f, MAE=%.3f, rho=%.3f, median_feats=%.0f, stable(>50%%)=%d\n",
              pw, R2_cv, MAE, ifelse(is.na(rho), 0, rho),
              median(fold_n_sig), length(stable_features)))
}

saveRDS(nested_cv_results, file.path(RESULTS_DIR, "nested_cv_results.rds"))

# Summary table
nested_summary <- data.frame(
  pathway = names(nested_cv_results),
  pathway_name = pw_lookup$pathway_name[match(names(nested_cv_results), pw_lookup$pw_col)],
  R2_cv = sapply(nested_cv_results, function(x) x$R2_cv),
  MAE = sapply(nested_cv_results, function(x) x$MAE),
  spearman = sapply(nested_cv_results, function(x) ifelse(is.na(x$spearman), 0, x$spearman)),
  n_stable_features = sapply(nested_cv_results, function(x) length(x$stable_features)),
  median_features_per_fold = sapply(nested_cv_results, function(x) median(x$fold_n_sig)),
  stringsAsFactors = FALSE
) %>% arrange(desc(R2_cv))

cat("\n--- Nested CV Summary (sorted by R2_cv) ---\n")
print(nested_summary, row.names = FALSE)

# Identify pathways with positive R2_cv
positive_pathways <- nested_summary$pathway[nested_summary$R2_cv > 0]
cat(sprintf("\nPathways with R2_cv > 0: %d (%s)\n",
            length(positive_pathways),
            ifelse(length(positive_pathways) > 0,
                   paste(positive_pathways, collapse = ", "), "none")))

write.csv(nested_summary, file.path(RESULTS_DIR, "nested_cv_summary.csv"),
          row.names = FALSE)

# ============================================================
# STEP 2: BOOTSTRAP CIs FOR NESTED CV
# ============================================================

cat("\n========== STEP 2: BOOTSTRAP CIs ==========\n")

B_BOOT <- 1000

for (pw in names(nested_cv_results)) {
  obs  <- nested_cv_results[[pw]]$observed
  pred <- nested_cv_results[[pw]]$predictions
  N    <- length(obs)

  boot_R2  <- numeric(B_BOOT)
  boot_MAE <- numeric(B_BOOT)

  set.seed(3000 + which(names(nested_cv_results) == pw))

  for (b in 1:B_BOOT) {
    idx <- sample(1:N, N, replace = TRUE)
    y_b <- obs[idx]
    p_b <- pred[idx]

    SS_res_b <- sum((y_b - p_b)^2)
    SS_tot_b <- sum((y_b - mean(obs))^2)

    boot_R2[b]  <- ifelse(SS_tot_b > 0, 1 - SS_res_b / SS_tot_b, NA)
    boot_MAE[b] <- mean(abs(y_b - p_b))
  }

  R2_ci  <- quantile(boot_R2, c(0.025, 0.975), na.rm = TRUE)
  MAE_ci <- quantile(boot_MAE, c(0.025, 0.975), na.rm = TRUE)

  nested_cv_results[[pw]]$R2_ci  <- R2_ci
  nested_cv_results[[pw]]$MAE_ci <- MAE_ci
}

# Report CIs for positive pathways
if (length(positive_pathways) > 0) {
  cat("\nBootstrap CIs for pathways with R2_cv > 0:\n")
  for (pw in positive_pathways) {
    r <- nested_cv_results[[pw]]
    cat(sprintf("  %s: R2=%.3f [%.3f, %.3f], MAE=%.3f [%.3f, %.3f]\n",
                pw, r$R2_cv, r$R2_ci[1], r$R2_ci[2],
                r$MAE, r$MAE_ci[1], r$MAE_ci[2]))
  }
} else {
  cat("  No pathways with R2_cv > 0 — reporting CIs for top 3:\n")
  top3 <- head(nested_summary$pathway, 3)
  for (pw in top3) {
    r <- nested_cv_results[[pw]]
    cat(sprintf("  %s: R2=%.3f [%.3f, %.3f]\n",
                pw, r$R2_cv, r$R2_ci[1], r$R2_ci[2]))
  }
}

# Update saved results with CIs
saveRDS(nested_cv_results, file.path(RESULTS_DIR, "nested_cv_results.rds"))

# ============================================================
# STEP 3: NESTED PERMUTATION TEST
# ============================================================

cat("\n========== STEP 3: NESTED PERMUTATION TEST ==========\n")

N_PERM_NESTED <- 1000
nested_perm_results <- list()

if (length(positive_pathways) == 0) {
  cat("  No pathways with R2_cv > 0 — skipping permutation test.\n")
  cat("  All nested CV R2 values are <= 0, so permutation p-values are trivially 1.0.\n")
} else {
  for (pw in positive_pathways) {
    cat(sprintf("\n  Permutation test for %s (N=%d perms)...\n", pw, N_PERM_NESTED))

    set.seed(5000 + which(pathway_cols == pw))

    z_pathway <- paste0("z_", pw)
    y <- patient_level[[z_pathway]]
    X_all <- as.matrix(patient_level[, all_z_features, drop = FALSE])
    N <- length(y)

    observed_R2 <- nested_cv_results[[pw]]$R2_cv
    perm_R2 <- rep(NA_real_, N_PERM_NESTED)

    t_start <- Sys.time()

    for (perm in 1:N_PERM_NESTED) {
      y_perm <- sample(y)
      perm_preds <- numeric(N)

      for (i in 1:N) {
        X_train <- X_all[-i, , drop = FALSE]
        y_train <- y_perm[-i]
        X_test  <- X_all[i, , drop = FALSE]

        # Feature selection inside fold
        spearman_p <- apply(X_train, 2, function(col) {
          tryCatch(cor.test(col, y_train, method = "spearman", exact = FALSE)$p.value,
                   error = function(e) 1)
        })
        fdr_vals <- p.adjust(spearman_p, method = "BH")
        sig_feats <- names(fdr_vals)[fdr_vals < 0.10]

        if (length(sig_feats) > 5) {
          sig_feats <- names(sort(spearman_p[sig_feats]))[1:5]
        }

        if (length(sig_feats) == 0) {
          perm_preds[i] <- mean(y_train)
          next
        }

        X_sel_train <- X_train[, sig_feats, drop = FALSE]
        X_sel_test  <- X_test[, sig_feats, drop = FALSE]

        best_alpha <- 0.5; best_lambda <- NULL; best_cvm <- Inf
        for (av in seq(0.1, 1.0, by = 0.1)) {
          cvf <- tryCatch(
            suppressWarnings(cv.glmnet(X_sel_train, y_train, alpha = av,
                                       nfolds = min(5, N - 1))),
            error = function(e) NULL)
          if (!is.null(cvf)) {
            idx_1se <- which(cvf$lambda == cvf$lambda.1se)
            if (length(idx_1se) > 0 && cvf$cvm[idx_1se] < best_cvm) {
              best_cvm <- cvf$cvm[idx_1se]; best_alpha <- av
              best_lambda <- cvf$lambda.1se
            }
          }
        }

        if (!is.null(best_lambda)) {
          fit <- glmnet(X_sel_train, y_train, alpha = best_alpha, lambda = best_lambda)
          perm_preds[i] <- predict(fit, newx = X_sel_test)[1]
        } else {
          perm_preds[i] <- mean(y_train)
        }
      }

      SS_res_p <- sum((y_perm - perm_preds)^2)
      SS_tot_p <- sum((y_perm - mean(y_perm))^2)
      perm_R2[perm] <- ifelse(SS_tot_p > 0, 1 - SS_res_p / SS_tot_p, NA)

      if (perm %% 50 == 0) {
        elapsed <- difftime(Sys.time(), t_start, units = "mins")
        rate <- perm / as.numeric(elapsed)
        remaining <- (N_PERM_NESTED - perm) / rate
        cat(sprintf("    perm %d/%d (%.1f min elapsed, ~%.1f min remaining)\n",
                    perm, N_PERM_NESTED, elapsed, remaining))
      }
    }

    perm_p <- mean(perm_R2 >= observed_R2, na.rm = TRUE)

    nested_perm_results[[pw]] <- list(
      observed_R2 = observed_R2,
      perm_R2     = perm_R2,
      perm_p      = perm_p,
      n_valid     = sum(!is.na(perm_R2))
    )

    cat(sprintf("  %s: R2_cv=%.3f, nested_perm_p=%.4f (%d valid perms, %.1f min)\n",
                pw, observed_R2, perm_p, sum(!is.na(perm_R2)),
                difftime(Sys.time(), t_start, units = "mins")))
  }
}

saveRDS(nested_perm_results, file.path(RESULTS_DIR, "nested_permutation_results.rds"))

# ============================================================
# STEP 4: CLINICAL COVARIATE ADJUSTMENT (S9)
# ============================================================

cat("\n========== STEP 4: CLINICAL COVARIATE ADJUSTMENT (S9) ==========\n")

tumor_details <- readRDS(file.path(DATA_DIR, "raw", "rnaseq", "ivygap_tumor_details.rds"))

# Create patient_id from tumor_name (W10-1-1 -> W10)
tumor_details$patient_id <- sub("-.*", "", tumor_details$tumor_name)

# Parse age
tumor_details$age <- as.numeric(sub(" yrs", "", tumor_details$age_in_years))

# Deduplicate (keep first entry per patient)
tumor_details <- tumor_details[!duplicated(tumor_details$patient_id), ]
cat(sprintf("Tumor details: %d unique patients\n", nrow(tumor_details)))

# Merge with analysis_data (observation-level for LMM)
analysis_data_s9 <- merge(
  analysis_data,
  tumor_details[, c("patient_id", "age", "mgmt_methylation", "molecular_subtype")],
  by = "patient_id", all.x = TRUE
)

cat(sprintf("After merge: %d rows, %d patients with covariates\n",
            nrow(analysis_data_s9),
            sum(!is.na(analysis_data_s9$age))))

cat(sprintf("  Age: %d non-missing (range: %.0f-%.0f)\n",
            sum(!is.na(analysis_data_s9$age)),
            min(analysis_data_s9$age, na.rm = TRUE),
            max(analysis_data_s9$age, na.rm = TRUE)))
cat(sprintf("  MGMT: %s\n",
            paste(names(table(analysis_data_s9$mgmt_methylation)),
                  table(analysis_data_s9$mgmt_methylation),
                  sep = "=", collapse = ", ")))
cat(sprintf("  Molecular subtype: %s\n",
            paste(names(table(analysis_data_s9$molecular_subtype)),
                  table(analysis_data_s9$molecular_subtype),
                  sep = "=", collapse = ", ")))

# Simplify molecular subtype (comma-separated -> Mixed)
analysis_data_s9$mol_subtype_simple <- ifelse(
  grepl(",", analysis_data_s9$molecular_subtype),
  "Mixed",
  analysis_data_s9$molecular_subtype
)

# Standardize age
analysis_data_s9$z_age <- as.numeric(scale(analysis_data_s9$age))

# Set subcompartment factor
analysis_data_s9$subcompartment <- factor(
  analysis_data_s9$subcompartment,
  levels = c("ET", "NET", "ED")
)

# Run S9 for pw_4 (Inflammatory Response)
s9_results <- list()

if ("pw_4" %in% names(lmm_results)) {
  top5 <- lmm_results[["pw_4"]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- "z_pw_4"

  # Define progressive covariate models
  models_s9 <- list(
    "A: Primary (no covariates)" = list(
      full = as.formula(paste(z_pathway, "~",
        paste(c(z_features, "subcompartment"), collapse = " + "),
        "+ (1 | patient_id)")),
      null = as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))
    ),
    "B: + age" = list(
      full = as.formula(paste(z_pathway, "~",
        paste(c(z_features, "subcompartment", "z_age"), collapse = " + "),
        "+ (1 | patient_id)")),
      null = as.formula(paste(z_pathway, "~ subcompartment + z_age + (1 | patient_id)"))
    ),
    "C: + age + MGMT" = list(
      full = as.formula(paste(z_pathway, "~",
        paste(c(z_features, "subcompartment", "z_age", "mgmt_methylation"), collapse = " + "),
        "+ (1 | patient_id)")),
      null = as.formula(paste(z_pathway,
        "~ subcompartment + z_age + mgmt_methylation + (1 | patient_id)"))
    ),
    "D: + age + MGMT + subtype" = list(
      full = as.formula(paste(z_pathway, "~",
        paste(c(z_features, "subcompartment", "z_age", "mgmt_methylation",
                "mol_subtype_simple"), collapse = " + "),
        "+ (1 | patient_id)")),
      null = as.formula(paste(z_pathway,
        "~ subcompartment + z_age + mgmt_methylation + mol_subtype_simple + (1 | patient_id)"))
    )
  )

  # Filter to complete cases
  covar_cols <- c(z_pathway, z_features, "subcompartment", "patient_id",
                  "z_age", "mgmt_methylation", "mol_subtype_simple")
  data_cc <- analysis_data_s9[complete.cases(analysis_data_s9[, covar_cols]), ]
  cat(sprintf("\nComplete cases for S9: %d rows, %d patients\n",
              nrow(data_cc), length(unique(data_cc$patient_id))))

  for (mod_name in names(models_s9)) {
    tryCatch({
      suppressWarnings(suppressMessages({
        full_ml <- lmer(models_s9[[mod_name]]$full, data = data_cc, REML = FALSE)
        null_ml <- lmer(models_s9[[mod_name]]$null, data = data_cc, REML = FALSE)
        lrt <- anova(null_ml, full_ml)

        full_reml <- lmer(models_s9[[mod_name]]$full, data = data_cc, REML = TRUE)
        r2 <- r2_nakagawa(full_reml, verbose = FALSE)
        coefs <- tidy(full_reml, effects = "fixed", conf.int = TRUE)
      }))

      # Extract radiomic feature coefficients only
      rad_coefs <- coefs %>%
        filter(grepl("^z_rad_", term)) %>%
        mutate(p_holm = p.adjust(p.value, method = "holm"))

      s9_results[[mod_name]] <- list(
        R2_marginal    = r2$R2_marginal,
        R2_conditional = r2$R2_conditional,
        LRT_p          = lrt$`Pr(>Chisq)`[2],
        LRT_chisq      = lrt$Chisq[2],
        coefficients   = coefs,
        rad_coefs      = rad_coefs,
        N_obs          = nrow(data_cc),
        N_patients     = length(unique(data_cc$patient_id)),
        converged      = !isSingular(full_reml)
      )

      cat(sprintf("  %s: R2_m=%.3f, R2_c=%.3f, LRT_p=%.4f, N=%d patients\n",
                  mod_name, r2$R2_marginal, r2$R2_conditional,
                  lrt$`Pr(>Chisq)`[2],
                  length(unique(data_cc$patient_id))))
    }, error = function(e) {
      cat(sprintf("  %s: FAILED - %s\n", mod_name, e$message))
      s9_results[[mod_name]] <<- list(
        R2_marginal = NA, LRT_p = NA,
        error = e$message
      )
    })
  }
}

saveRDS(s9_results, file.path(RESULTS_DIR, "clinical_covariates_s9.rds"))

# ============================================================
# STEP 5: EXPORT FEATURE LOOKUP TABLE
# ============================================================

cat("\n========== STEP 5: EXPORT FEATURE LOOKUP ==========\n")

# Collect all features used across nested CV and LMM
all_nested_features <- unique(unlist(lapply(nested_cv_results, function(x) {
  unique(unlist(x$fold_features))
})))
all_nested_features_raw <- sub("^z_", "", all_nested_features)

lmm_features_raw <- unique(unlist(lapply(lmm_results, function(x) x$feature_names)))
all_used_features <- unique(c(all_nested_features_raw, lmm_features_raw))

# Build enriched lookup table
feature_lookup_table <- rad_lookup[rad_lookup$rad_col %in% all_used_features, ]
feature_lookup_table <- feature_lookup_table %>%
  mutate(
    sequence = sub("_.*", "", feature_name),
    feature_type = case_when(
      grepl("CoLlAGe|CoLIAGe", feature_name) ~ "Higher-order (CoLIAGe)",
      grepl("GLSZM", feature_name)            ~ "Texture (GLSZM)",
      grepl("GLCM", feature_name)             ~ "Texture (GLCM)",
      grepl("GLRLM", feature_name)            ~ "Texture (GLRLM)",
      grepl("Histogram", feature_name)        ~ "First-order (Histogram)",
      grepl("Laws", feature_name)             ~ "Texture (Laws)",
      grepl("Gabor", feature_name)            ~ "Texture (Gabor)",
      grepl("Haralick", feature_name)         ~ "Texture (Haralick)",
      grepl("NGTDM", feature_name)            ~ "Texture (NGTDM)",
      grepl("shape|Shape", feature_name)      ~ "Shape",
      TRUE                                    ~ "Other"
    )
  ) %>%
  arrange(rad_col)

write.csv(feature_lookup_table, file.path(RESULTS_DIR, "feature_lookup.csv"),
          row.names = FALSE)
cat(sprintf("Feature lookup exported: %d features used across all analyses\n",
            nrow(feature_lookup_table)))

# Export pw_4 stability results with IBSI names
if ("pw_4" %in% names(nested_cv_results)) {
  pw4_stability <- nested_cv_results[["pw_4"]]$feat_stability
  if (nrow(pw4_stability) > 0) {
    pw4_stability$rad_col <- sub("^z_", "", pw4_stability$feature)
    pw4_stability <- merge(pw4_stability, rad_lookup, by = "rad_col", all.x = TRUE)
    pw4_stability <- pw4_stability %>% arrange(desc(n_folds))
    write.csv(pw4_stability, file.path(RESULTS_DIR, "pw4_feature_stability.csv"),
              row.names = FALSE)
    cat(sprintf("pw_4 feature stability: %d features ever selected, %d stable (>50%%)\n",
                nrow(pw4_stability),
                sum(pw4_stability$pct_folds > 50)))

    # Print top features for pw_4
    cat("\nTop features for pw_4 (Inflammatory Response) by stability:\n")
    top_feats <- head(pw4_stability, 10)
    for (j in seq_len(nrow(top_feats))) {
      cat(sprintf("  %s (%s): %d/%d folds (%.0f%%)\n",
                  top_feats$rad_col[j],
                  top_feats$feature_name[j],
                  top_feats$n_folds[j], N_patients,
                  top_feats$pct_folds[j]))
    }
  }
}

# ============================================================
# STEP 6: COMPARISON TABLE (NESTED CV vs LEGACY)
# ============================================================

cat("\n========== STEP 6: NESTED CV vs LEGACY COMPARISON ==========\n")

comparison <- data.frame(
  pathway      = character(),
  pathway_name = character(),
  legacy_prescreen_R2 = numeric(),
  legacy_s5_R2        = numeric(),
  nested_R2_cv        = numeric(),
  stringsAsFactors = FALSE
)

for (pw in names(nested_cv_results)) {
  # Legacy pre-screened EN result (from original pipeline)
  legacy_prescreen <- if (pw %in% names(elastic_net_results)) {
    elastic_net_results[[pw]]$R2_cv
  } else {
    NA
  }

  # Legacy S5 result (feature selection inside CV, but only for target_pathways)
  legacy_s5 <- if (pw %in% names(en_results_s5)) {
    en_results_s5[[pw]]$R2_cv
  } else {
    NA
  }

  nested_R2 <- nested_cv_results[[pw]]$R2_cv

  comparison <- rbind(comparison, data.frame(
    pathway             = pw,
    pathway_name        = pw_lookup$pathway_name[pw_lookup$pw_col == pw],
    legacy_prescreen_R2 = legacy_prescreen,
    legacy_s5_R2        = legacy_s5,
    nested_R2_cv        = nested_R2,
    stringsAsFactors = FALSE
  ))
}

comparison <- comparison %>% arrange(desc(nested_R2_cv))

cat("\n")
print(comparison, row.names = FALSE)

write.csv(comparison, file.path(RESULTS_DIR, "nested_vs_legacy_comparison.csv"),
          row.names = FALSE)

# ============================================================
# STEP 7: FINAL SUMMARY
# ============================================================

cat("\n========== FINAL SUMMARY ==========\n")

cat("\n1. NESTED CV RESULTS (all 24 pathways):\n")
cat(sprintf("   Pathways with R2_cv > 0: %d\n", length(positive_pathways)))
if (length(positive_pathways) > 0) {
  for (pw in positive_pathways) {
    r <- nested_cv_results[[pw]]
    cat(sprintf("   - %s (%s): R2_cv=%.3f [%.3f, %.3f]\n",
                pw, pw_lookup$pathway_name[pw_lookup$pw_col == pw],
                r$R2_cv, r$R2_ci[1], r$R2_ci[2]))
  }
}

cat("\n2. NESTED PERMUTATION RESULTS:\n")
if (length(nested_perm_results) > 0) {
  for (pw in names(nested_perm_results)) {
    r <- nested_perm_results[[pw]]
    cat(sprintf("   - %s: R2_cv=%.3f, perm_p=%.4f\n", pw, r$observed_R2, r$perm_p))
  }
} else {
  cat("   No pathways with R2_cv > 0 — all permutation p-values = 1.0\n")
}

cat("\n3. CLINICAL COVARIATE ADJUSTMENT (S9):\n")
if (length(s9_results) > 0) {
  for (mod_name in names(s9_results)) {
    r <- s9_results[[mod_name]]
    if (!is.na(r$R2_marginal)) {
      cat(sprintf("   %s: R2_m=%.3f, LRT_p=%.4f\n", mod_name, r$R2_marginal, r$LRT_p))
    }
  }
}

cat("\n4. KEY COMPARISON:\n")
if ("pw_4" %in% names(nested_cv_results)) {
  legacy_pre <- if ("pw_4" %in% names(elastic_net_results)) {
    elastic_net_results[["pw_4"]]$R2_cv
  } else NA
  legacy_s5v <- if ("pw_4" %in% names(en_results_s5)) {
    en_results_s5[["pw_4"]]$R2_cv
  } else NA
  nested_v <- nested_cv_results[["pw_4"]]$R2_cv

  cat(sprintf("   pw_4 (Inflammatory Response):\n"))
  cat(sprintf("     Legacy pre-screened EN: R2_cv = %.3f\n",
              ifelse(is.na(legacy_pre), NA, legacy_pre)))
  cat(sprintf("     Legacy S5 (inside CV):  R2_cv = %.3f\n",
              ifelse(is.na(legacy_s5v), NA, legacy_s5v)))
  cat(sprintf("     Nested CV (this script): R2_cv = %.3f\n", nested_v))
}

# ============================================================
# SAVE SESSION
# ============================================================

cat("\n========== SAVING SESSION ==========\n")

save.image(file.path(RESULTS_DIR, "03b_nested_cv_session.RData"))

sink(file.path(RESULTS_DIR, "03b_session_info.txt"))
cat("Session info for 03b_nested_cv_analysis.R\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Seed: 42\n\n"))
sessionInfo()
sink()

cat(sprintf("Session saved to: %s\n",
            file.path(RESULTS_DIR, "03b_nested_cv_session.RData")))
cat(sprintf("\n=== SCRIPT 03b COMPLETE (%s) ===\n", Sys.time()))
