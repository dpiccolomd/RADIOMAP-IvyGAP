# ============================================================
# RADIOMAP-IvyGAP: Statistical Analysis Pipeline
# Script 03: Mixed-Effects Models + Elastic Net + Permutation
# ============================================================
# Author: Daniele Piccolo, MD
# Date: 2026-02-21
# Dependencies: lme4, lmerTest, performance, glmnet, pmsampsize,
#               GSVA, msigdbr, caret, dplyr, tidyr, broom.mixed,
#               ggplot2, patchwork, pbkrtest (optional, for KR df)
# Input: Prepared data from scripts 01 (data audit) and 02 (ssGSEA)
# Output: Association results, prediction results, figures, tables
# ============================================================

# --- Setup ---
library(lme4)
library(lmerTest)
library(performance)
library(glmnet)
library(pmsampsize)
library(caret)
library(dplyr)
library(tidyr)
library(broom.mixed)
library(ggplot2)
library(patchwork)

set.seed(42)

PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# STEP 0: POWER ANALYSIS (pmsampsize)
# ============================================================

cat("\n========== STEP 0: POWER ANALYSIS ==========\n")

# Riley et al. 2020 — continuous outcome, criterion 1 (shrinkage)
pms_result <- pmsampsize(
  type       = "c",
  rsquared   = 0.20,
  parameters = 5,
  shrinkage  = 0.90,
  intercept  = 0,
  sd         = 1
)
cat("Riley 2020 result for P=5, R2=0.20, S>=0.90:\n")
print(pms_result)

# Minimum detectable R2 at available N
for (N in c(25, 28, 31)) {
  # From Riley C1: S = 1 - (P/N)*((1-R2)/R2)  [van Houwelingen shrinkage]
  # Rearranging for R2 at fixed N, S=0.90, P=5:
  #   (1-S) = (P/N)*((1-R2)/R2)
  #   R2/(1-R2) = P / (N*(1-S))
  #   R2 = P / (N*(1-S) + P)
  R2_min <- 5 / (N * (1 - 0.90) + 5)
  cat(sprintf("N=%d: Min detectable R2 for S>=0.90 = %.3f\n", N, R2_min))
}

# --- Attenuation factor (Park 2025 mapping accuracy) ---
cat("\n--- Attenuation Factor (Park 2025) ---\n")
park_r <- c(CT_ET = 0.238, PAN_NET = 0.241, MVP_ET = 0.195, IT_LE_ED = 0.294)
mean_r <- mean(park_r)
attenuation <- mean_r^2
cat(sprintf("Park 2025 zone-subcompartment correlations: %s\n",
            paste(sprintf("%s=%.3f", names(park_r), park_r), collapse = ", ")))
cat(sprintf("Mean r = %.3f, attenuation factor (r^2) = %.3f\n", mean_r, attenuation))
cat(sprintf("Implication: observed associations are attenuated by ~%.0f%% due to spatial mismatch\n",
            (1 - attenuation) * 100))
cat(sprintf("To detect rho_true = 0.50, need observed r >= %.3f at this attenuation\n",
            0.50 * mean_r))

# ============================================================
# STEP 1: LOAD PREPARED DATA
# ============================================================
# Expected input: analysis_data from script 02
# Columns: patient_id, subcompartment, pathway scores, radiomic features
# One row per patient x subcompartment (3 rows per patient)

cat("\n========== STEP 1: LOAD DATA ==========\n")

# Load analysis-ready dataset (output from script 02_ssgsea_aggregation.R)
analysis_data <- readRDS(file.path(DATA_DIR, "analysis_data.rds"))

# Expected structure:
# - patient_id: character (e.g., "W1", "W2", ...)
# - subcompartment: factor(ET, NET, ED)
# - pathway columns: ssGSEA scores (24 pathways)
# - radiomic columns: feature values (reduced from 11,700)

cat(sprintf("Loaded: %d rows, %d columns\n", nrow(analysis_data), ncol(analysis_data)))
cat(sprintf("Patients: %d\n", length(unique(analysis_data$patient_id))))
cat(sprintf("Subcompartments per patient: %s\n",
            paste(table(table(analysis_data$patient_id)), collapse = ", ")))

# Define column groups (to be updated with actual column names)
pathway_cols  <- grep("^pw_", colnames(analysis_data), value = TRUE)
radiomic_cols <- grep("^rad_", colnames(analysis_data), value = TRUE)

cat(sprintf("Pathway columns: %d\n", length(pathway_cols)))
cat(sprintf("Radiomic columns: %d\n", length(radiomic_cols)))

# Validate expected structure
stopifnot(
  "patient_id column missing" = "patient_id" %in% colnames(analysis_data),
  "subcompartment column missing" = "subcompartment" %in% colnames(analysis_data),
  "Expected 24 pathway columns" = length(pathway_cols) == 24,
  "No radiomic columns found" = length(radiomic_cols) > 0,
  "Expected >= 25 patients" = length(unique(analysis_data$patient_id)) >= 25
)

# ============================================================
# STEP 2: FEATURE REDUCTION
# ============================================================

cat("\n========== STEP 2: FEATURE REDUCTION ==========\n")

# 2a. Near-zero variance filter (caret::nearZeroVar)
# Removes features with very low variance relative to their own distribution,
# using frequency ratio and percent-unique criteria (scale-invariant)
nzv_idx <- nearZeroVar(analysis_data[, radiomic_cols])
if (length(nzv_idx) > 0) {
  features_pass_nzv <- radiomic_cols[-nzv_idx]
} else {
  features_pass_nzv <- radiomic_cols
}
cat(sprintf("After near-zero-variance filter: %d / %d features remain\n",
            length(features_pass_nzv), length(radiomic_cols)))

# 2b. High correlation filter (|Spearman r| > 0.90)
# Plan specifies Spearman correlation for clustering; caret::findCorrelation
# accepts any pre-computed correlation matrix
cor_matrix <- cor(analysis_data[, features_pass_nzv], use = "complete.obs",
                  method = "spearman")
high_cor_idx <- findCorrelation(cor_matrix, cutoff = 0.90, names = TRUE)
features_pass_cor <- setdiff(features_pass_nzv, high_cor_idx)
cat(sprintf("After correlation filter (|r|>0.90): %d features remain\n",
            length(features_pass_cor)))

# 2c. Univariate Spearman filter per pathway — PER SUBCOMPARTMENT (plan spec)
# Run correlations within each subcompartment, then union significant features.
# Top 5 selected by minimum FDR across subcompartments.
univariate_results <- list()

subcompartments <- unique(analysis_data$subcompartment)

for (pw in pathway_cols) {
  # Collect per-subcompartment results
  all_sc_results <- list()

  for (sc in subcompartments) {
    sc_data <- analysis_data[analysis_data$subcompartment == sc, ]

    spearman_res <- data.frame(
      feature        = features_pass_cor,
      subcompartment = sc,
      rho            = numeric(length(features_pass_cor)),
      p_raw          = numeric(length(features_pass_cor)),
      stringsAsFactors = FALSE
    )

    for (i in seq_along(features_pass_cor)) {
      feat <- features_pass_cor[i]
      test <- cor.test(
        sc_data[[feat]],
        sc_data[[pw]],
        method = "spearman",
        exact  = FALSE
      )
      spearman_res$rho[i]   <- test$estimate
      spearman_res$p_raw[i] <- test$p.value
    }

    spearman_res$p_fdr <- p.adjust(spearman_res$p_raw, method = "BH")
    all_sc_results[[sc]] <- spearman_res
  }

  # Combine across subcompartments
  combined <- do.call(rbind, all_sc_results)

  # For each feature, take its best (minimum) FDR across subcompartments
  best_per_feature <- combined %>%
    group_by(feature) %>%
    summarise(
      min_fdr  = min(p_fdr, na.rm = TRUE),
      best_rho = rho[which.min(p_fdr)],
      best_sc  = subcompartment[which.min(p_fdr)],
      .groups  = "drop"
    ) %>%
    filter(min_fdr < 0.10) %>%
    arrange(min_fdr)

  top5 <- head(best_per_feature$feature, 5)

  univariate_results[[pw]] <- list(
    per_subcompartment = all_sc_results,
    best_per_feature   = best_per_feature,
    top_5              = top5,
    n_significant      = nrow(best_per_feature)
  )

  cat(sprintf("  %s: %d features pass FDR<0.10 (any SC), top %d selected\n",
              pw, nrow(best_per_feature), length(top5)))
}

# Save feature reduction results
saveRDS(univariate_results, file.path(RESULTS_DIR, "univariate_results.rds"))

# ============================================================
# STEP 3: STANDARDIZE FEATURES
# ============================================================

cat("\n========== STEP 3: STANDARDIZE ==========\n")

# Z-score radiomic features WITHIN each subcompartment
analysis_data <- analysis_data %>%
  group_by(subcompartment) %>%
  mutate(across(all_of(features_pass_cor),
                ~ as.numeric(scale(.)),
                .names = "z_{.col}")) %>%
  ungroup()

# Z-score pathway scores globally
analysis_data <- analysis_data %>%
  mutate(across(all_of(pathway_cols),
                ~ as.numeric(scale(.)),
                .names = "z_{.col}"))

# Set subcompartment factor levels (ET = reference)
analysis_data$subcompartment <- factor(
  analysis_data$subcompartment,
  levels = c("ET", "NET", "ED")
)

cat("Standardization complete.\n")

# ============================================================
# STEP 4: PRIMARY ANALYSIS — MIXED-EFFECTS MODELS
# ============================================================

cat("\n========== STEP 4: MIXED-EFFECTS MODELS ==========\n")

# Helper: check lmer convergence
check_convergence <- function(model) {
  cc <- model@optinfo$conv$lme4
  has_warnings <- !is.null(cc$messages) && length(cc$messages) > 0
  is_singular <- isSingular(model)
  list(converged = !has_warnings && !is_singular,
       singular  = is_singular,
       messages  = cc$messages)
}

lmm_results <- list()
model_p_values <- numeric()

for (pw in pathway_cols) {
  top5 <- univariate_results[[pw]]$top_5

  if (length(top5) < 1) {
    cat(sprintf("  SKIP %s: no features pass univariate filter\n", pw))
    next
  }

  n_feat <- min(length(top5), 5)
  z_features <- paste0("z_", top5[1:n_feat])
  z_pathway  <- paste0("z_", pw)

  # --- Full model (Option B) ---
  fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
  full_formula <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))

  # --- Null model (subcompartment only) ---
  null_formula <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

  # Fit both with ML (required for LRT)
  full_model_ml <- lmer(full_formula, data = analysis_data, REML = FALSE)
  null_model_ml <- lmer(null_formula, data = analysis_data, REML = FALSE)

  # Likelihood Ratio Test
  lrt <- anova(null_model_ml, full_model_ml)
  lrt_p <- lrt$`Pr(>Chisq)`[2]

  # Re-fit full model with REML for final estimates
  full_model_reml <- lmer(full_formula, data = analysis_data, REML = TRUE)

  # R2 decomposition
  r2_vals <- r2_nakagawa(full_model_reml, verbose = FALSE)

  # Coefficient table
  coef_tbl <- tidy(full_model_reml, effects = "fixed", conf.int = TRUE)

  # ANOVA Type II (appropriate with treatment coding; Type III requires sum contrasts)
  anova_tbl <- anova(full_model_reml, type = 2, ddf = "Satterthwaite")

  # Check convergence
  conv_info <- check_convergence(full_model_reml)
  if (!conv_info$converged) {
    cat(sprintf("  WARNING %s: convergence issue — singular=%s, messages=%s\n",
                pw, conv_info$singular,
                paste(conv_info$messages, collapse = "; ")))
  }

  # Store
  lmm_results[[pw]] <- list(
    model_reml     = full_model_reml,
    formula_str    = deparse(full_formula),
    n_features     = n_feat,
    feature_names  = top5[1:n_feat],
    R2_marginal    = r2_vals$R2_marginal,
    R2_conditional = r2_vals$R2_conditional,
    LRT_p          = lrt_p,
    anova          = anova_tbl,
    coefficients   = coef_tbl,
    AIC            = AIC(full_model_reml),
    BIC            = BIC(full_model_reml),
    converged      = conv_info$converged,
    singular       = conv_info$singular
  )

  model_p_values[pw] <- lrt_p

  cat(sprintf("  %s: R2_m=%.3f, R2_c=%.3f, LRT_p=%.4f, P=%d\n",
              pw, r2_vals$R2_marginal, r2_vals$R2_conditional, lrt_p, n_feat))
}

# --- Level 1: FDR correction across pathways ---
model_fdr <- p.adjust(model_p_values, method = "BH")

results_summary <- data.frame(
  pathway        = names(lmm_results),
  R2_marginal    = sapply(lmm_results, function(x) x$R2_marginal),
  R2_conditional = sapply(lmm_results, function(x) x$R2_conditional),
  p_LRT          = model_p_values[names(lmm_results)],
  p_FDR          = model_fdr[names(lmm_results)],
  n_features     = sapply(lmm_results, function(x) x$n_features),
  stringsAsFactors = FALSE
) %>%
  arrange(p_FDR)

cat("\n--- Results Summary (sorted by FDR) ---\n")
print(results_summary, row.names = FALSE)

# Identify significant pathways
sig_pathways <- results_summary %>% filter(p_FDR < 0.05)
cat(sprintf("\nSignificant pathways (FDR<0.05): %d / %d\n",
            nrow(sig_pathways), nrow(results_summary)))

# --- Level 2: Coefficient-level tests within significant models ---
if (nrow(sig_pathways) > 0) {
  cat("\n--- Level 2: Coefficient-level results ---\n")
  for (pw in sig_pathways$pathway) {
    coefs <- lmm_results[[pw]]$coefficients %>%
      filter(!grepl("Intercept|subcompartment", term)) %>%
      mutate(p_holm = p.adjust(p.value, method = "holm"))

    cat(sprintf("\n  %s (FDR=%.4f, R2_m=%.3f):\n", pw,
                results_summary$p_FDR[results_summary$pathway == pw],
                lmm_results[[pw]]$R2_marginal))
    print(coefs[, c("term", "estimate", "conf.low", "conf.high", "p.value", "p_holm")])
  }
}

# Save results
saveRDS(lmm_results, file.path(RESULTS_DIR, "lmm_results.rds"))
write.csv(results_summary, file.path(RESULTS_DIR, "lmm_results_summary.csv"),
          row.names = FALSE)

# ============================================================
# STEP 5: PERMUTATION TEST FOR MIXED-EFFECTS MODELS
# ============================================================

cat("\n========== STEP 5: PERMUTATION TESTING (LMM) ==========\n")

N_PERM <- 1000
permutation_results <- list()

for (pw in names(lmm_results)) {
  # Per-pathway seed for reproducibility
  pw_idx <- which(names(lmm_results) == pw)
  set.seed(42 + pw_idx)

  observed_R2m <- lmm_results[[pw]]$R2_marginal

  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)
  fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
  full_formula <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))
  null_formula <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

  # Observed LRT chi-squared (test statistic: directly tests radiomic contribution)
  obs_full_ml <- lmer(full_formula, data = analysis_data, REML = FALSE)
  obs_null_ml <- lmer(null_formula, data = analysis_data, REML = FALSE)
  obs_lrt <- anova(obs_null_ml, obs_full_ml)
  observed_chisq <- obs_lrt$Chisq[2]

  perm_chisq <- rep(NA_real_, N_PERM)
  perm_R2m   <- rep(NA_real_, N_PERM)
  patient_ids <- unique(analysis_data$patient_id)

  for (i in 1:N_PERM) {
    # Permute at PATIENT level: shuffle which patient's pathway scores
    # get paired with which patient's radiomic features.
    # DESIGN NOTE: The random intercept (1|patient_id) retains the ORIGINAL
    # patient_id (radiomic data structure). This is intentional — we test whether
    # radiomic-pathway associations exist beyond patient-level clustering of
    # radiomic features. The pathway scores are shuffled between patients while
    # the radiomic grouping structure is preserved.
    shuffled_ids <- sample(patient_ids)
    id_map <- setNames(shuffled_ids, patient_ids)

    perm_data <- analysis_data
    perm_data$patient_id_perm <- id_map[perm_data$patient_id]

    # Get pathway scores indexed by permuted patient ID
    pathway_lookup <- analysis_data %>%
      select(patient_id, subcompartment, all_of(z_pathway)) %>%
      rename(patient_id_perm = patient_id)

    perm_data <- perm_data %>%
      select(-all_of(z_pathway)) %>%
      left_join(pathway_lookup, by = c("patient_id_perm", "subcompartment"))

    tryCatch({
      suppressWarnings(suppressMessages({
        # Use ML for LRT (primary test statistic)
        perm_full_ml <- lmer(full_formula, data = perm_data, REML = FALSE)
        perm_null_ml <- lmer(null_formula, data = perm_data, REML = FALSE)
        perm_lrt <- anova(perm_null_ml, perm_full_ml)
        perm_chisq[i] <- perm_lrt$Chisq[2]

        # Also store marginal R2 for secondary comparison
        perm_model_reml <- lmer(full_formula, data = perm_data, REML = TRUE)
        perm_R2m[i] <- r2_nakagawa(perm_model_reml, verbose = FALSE)$R2_marginal
      }))
    }, error = function(e) {
      perm_chisq[i] <<- NA
      perm_R2m[i]   <<- NA
    })
  }

  # Primary: LRT chi-squared permutation p-value
  perm_p_chisq <- mean(perm_chisq >= observed_chisq, na.rm = TRUE)
  # Secondary: marginal R2 permutation p-value
  perm_p_R2m <- mean(perm_R2m >= observed_R2m, na.rm = TRUE)

  permutation_results[[pw]] <- list(
    observed_R2m    = observed_R2m,
    observed_chisq  = observed_chisq,
    perm_chisq      = perm_chisq,
    perm_R2m        = perm_R2m,
    perm_p_chisq    = perm_p_chisq,
    perm_p_R2m      = perm_p_R2m,
    n_valid         = sum(!is.na(perm_chisq))
  )

  cat(sprintf("  %s: R2_m=%.3f, LRT_chisq=%.2f, perm_p(chisq)=%.4f, perm_p(R2m)=%.4f (%d valid)\n",
              pw, observed_R2m, observed_chisq, perm_p_chisq, perm_p_R2m,
              sum(!is.na(perm_chisq))))
}

saveRDS(permutation_results, file.path(RESULTS_DIR, "permutation_results_lmm.rds"))

# ============================================================
# STEP 6: SECONDARY ANALYSIS — ELASTIC NET WITH LOPO-CV
# ============================================================

cat("\n========== STEP 6: ELASTIC NET LOPO-CV ==========\n")

# Aggregate to patient level using RAW (non-z-scored) features and pathway scores.
# Z-scoring within-subcompartment values and then averaging across subcompartments
# is not meaningful — instead, average raw values then z-score at patient level.
raw_cols_to_avg <- c(pathway_cols, features_pass_cor)  # non z_ prefixed
patient_level <- analysis_data %>%
  group_by(patient_id) %>%
  summarise(across(all_of(raw_cols_to_avg), ~ mean(., na.rm = TRUE)),
            .groups = "drop")

# Z-score at patient level for Elastic Net
patient_level <- patient_level %>%
  mutate(across(all_of(features_pass_cor),
                ~ as.numeric(scale(.)),
                .names = "z_{.col}"))
patient_level <- patient_level %>%
  mutate(across(all_of(pathway_cols),
                ~ as.numeric(scale(.)),
                .names = "z_{.col}"))

N_patients <- nrow(patient_level)
cat(sprintf("Patient-level data: N=%d\n", N_patients))

# Run Elastic Net only for pathways that passed mixed-effects screen
# (or for all pathways if doing comprehensive analysis)
target_pathways <- if (nrow(sig_pathways) > 0) {
  sig_pathways$pathway
} else {
  # If no significant pathways, run for top 5 by raw p-value (exploratory)
  head(results_summary$pathway, 5)
}

elastic_net_results <- list()

for (pw in target_pathways) {
  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)

  y <- patient_level[[z_pathway]]
  X <- as.matrix(patient_level[, z_features, drop = FALSE])

  if (any(is.na(y)) || any(is.na(X))) {
    cat(sprintf("  SKIP %s: missing values\n", pw))
    next
  }

  N <- length(y)
  predictions <- numeric(N)
  best_alphas <- numeric(N)
  best_lambdas <- numeric(N)

  # ---- LOPO-CV ----
  for (i in 1:N) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X[i, , drop = FALSE]

    # Inner CV: grid search over alpha
    best_alpha  <- NA
    best_lambda <- NA
    best_cvm    <- Inf

    for (alpha_val in seq(0.1, 1.0, by = 0.1)) {
      cv_fit <- tryCatch({
        cv.glmnet(
          x            = X_train,
          y            = y_train,
          alpha        = alpha_val,
          nfolds       = min(5, N - 1),
          type.measure = "mse"
        )
      }, error = function(e) NULL)

      if (!is.null(cv_fit)) {
        # Use lambda.1se (more conservative, appropriate for small N)
        idx_1se <- which(cv_fit$lambda == cv_fit$lambda.1se)
        cvm_1se <- cv_fit$cvm[idx_1se]
        if (length(cvm_1se) > 0 && cvm_1se < best_cvm) {
          best_cvm    <- cvm_1se
          best_alpha  <- alpha_val
          best_lambda <- cv_fit$lambda.1se
        }
      }
    }

    # Predict left-out patient
    final_fit <- glmnet(X_train, y_train, alpha = best_alpha, lambda = best_lambda)
    predictions[i] <- predict(final_fit, newx = X_test)[1]
    best_alphas[i]  <- best_alpha
    best_lambdas[i] <- best_lambda
  }

  # ---- Performance metrics ----
  SS_res <- sum((y - predictions)^2)
  SS_tot <- sum((y - mean(y))^2)
  R2_cv  <- 1 - SS_res / SS_tot
  MAE    <- mean(abs(y - predictions))
  rho    <- cor(y, predictions, method = "spearman")

  # ---- Null model (intercept-only LOPO-CV) ----
  null_preds <- numeric(N)
  for (i in 1:N) {
    null_preds[i] <- mean(y[-i])
  }
  SS_res_null <- sum((y - null_preds)^2)
  R2_null <- 1 - SS_res_null / SS_tot

  elastic_net_results[[pw]] <- list(
    R2_cv       = R2_cv,
    MAE         = MAE,
    spearman    = rho,
    R2_null     = R2_null,
    predictions = predictions,
    observed    = y,
    best_alphas = best_alphas,
    best_lambdas = best_lambdas
  )

  cat(sprintf("  %s: R2_cv=%.3f, MAE=%.3f, rho=%.3f, R2_null=%.3f\n",
              pw, R2_cv, MAE, rho, R2_null))
}

# ---- Bootstrap CIs ----
cat("\n--- Bootstrap CIs for Elastic Net ---\n")
B_BOOT <- 1000

for (pw in names(elastic_net_results)) {
  obs  <- elastic_net_results[[pw]]$observed
  pred <- elastic_net_results[[pw]]$predictions
  N    <- length(obs)

  boot_R2  <- numeric(B_BOOT)
  boot_MAE <- numeric(B_BOOT)
  boot_rho <- numeric(B_BOOT)

  for (b in 1:B_BOOT) {
    idx <- sample(1:N, N, replace = TRUE)
    y_b <- obs[idx]
    p_b <- pred[idx]

    SS_res_b <- sum((y_b - p_b)^2)
    # Use original mean (not resampled mean) so R2 is comparable across bootstrap replicates
    SS_tot_b <- sum((y_b - mean(obs))^2)

    boot_R2[b]  <- ifelse(SS_tot_b > 0, 1 - SS_res_b / SS_tot_b, NA)
    boot_MAE[b] <- mean(abs(y_b - p_b))
    boot_rho[b] <- tryCatch(cor(y_b, p_b, method = "spearman"), error = function(e) NA)
  }

  R2_ci  <- quantile(boot_R2, c(0.025, 0.975), na.rm = TRUE)
  MAE_ci <- quantile(boot_MAE, c(0.025, 0.975), na.rm = TRUE)
  rho_ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)

  elastic_net_results[[pw]]$R2_ci  <- R2_ci
  elastic_net_results[[pw]]$MAE_ci <- MAE_ci
  elastic_net_results[[pw]]$rho_ci <- rho_ci

  cat(sprintf("  %s: R2=%.3f [%.3f, %.3f], MAE=%.3f [%.3f, %.3f], rho=%.3f [%.3f, %.3f]\n",
              pw,
              elastic_net_results[[pw]]$R2_cv, R2_ci[1], R2_ci[2],
              elastic_net_results[[pw]]$MAE, MAE_ci[1], MAE_ci[2],
              elastic_net_results[[pw]]$spearman, rho_ci[1], rho_ci[2]))
}

# ---- Optimism-corrected bootstrap (Harrell 2015) — PRIMARY CI method ----
cat("\n--- Optimism-corrected bootstrap for Elastic Net ---\n")
B_OPT <- 200

# Helper: fit EN with alpha grid, return model + apparent R2
fit_en_full <- function(X, y) {
  best_alpha <- 0.5; best_lambda <- NULL; best_cvm <- Inf
  for (av in seq(0.1, 1.0, by = 0.1)) {
    cvf <- tryCatch(cv.glmnet(X, y, alpha = av, nfolds = min(5, nrow(X))),
                    error = function(e) NULL)
    if (!is.null(cvf)) {
      idx <- which(cvf$lambda == cvf$lambda.1se)
      if (length(idx) > 0 && cvf$cvm[idx] < best_cvm) {
        best_cvm <- cvf$cvm[idx]; best_alpha <- av; best_lambda <- cvf$lambda.1se
      }
    }
  }
  if (is.null(best_lambda)) return(NULL)
  mod <- glmnet(X, y, alpha = best_alpha, lambda = best_lambda)
  list(model = mod, alpha = best_alpha, lambda = best_lambda)
}

for (pw in names(elastic_net_results)) {
  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)

  y <- patient_level[[z_pathway]]
  X <- as.matrix(patient_level[, z_features, drop = FALSE])
  N <- length(y)

  # Apparent R2: fit on all data, evaluate on all data
  full_fit <- fit_en_full(X, y)
  if (is.null(full_fit)) {
    cat(sprintf("  %s: SKIP optimism bootstrap (no model fit)\n", pw))
    next
  }
  apparent_pred <- predict(full_fit$model, newx = X)[, 1]
  apparent_R2 <- 1 - sum((y - apparent_pred)^2) / sum((y - mean(y))^2)

  set.seed(2000 + which(names(elastic_net_results) == pw))
  optimism <- rep(NA_real_, B_OPT)

  for (b in 1:B_OPT) {
    idx <- sample(1:N, N, replace = TRUE)
    X_boot <- X[idx, , drop = FALSE]
    y_boot <- y[idx]

    boot_fit <- tryCatch(fit_en_full(X_boot, y_boot), error = function(e) NULL)
    if (is.null(boot_fit)) next

    # Apparent performance on bootstrap sample
    p_boot <- predict(boot_fit$model, newx = X_boot)[, 1]
    ss_res_boot <- sum((y_boot - p_boot)^2)
    ss_tot_boot <- sum((y_boot - mean(y_boot))^2)
    R2_boot_app <- ifelse(ss_tot_boot > 0, 1 - ss_res_boot / ss_tot_boot, NA)

    # Test performance on original data
    p_orig <- predict(boot_fit$model, newx = X)[, 1]
    ss_res_orig <- sum((y - p_orig)^2)
    ss_tot_orig <- sum((y - mean(y))^2)
    R2_boot_test <- ifelse(ss_tot_orig > 0, 1 - ss_res_orig / ss_tot_orig, NA)

    optimism[b] <- R2_boot_app - R2_boot_test
  }

  corrected_R2 <- apparent_R2 - mean(optimism, na.rm = TRUE)

  elastic_net_results[[pw]]$R2_apparent          <- apparent_R2
  elastic_net_results[[pw]]$R2_optimism_corrected <- corrected_R2
  elastic_net_results[[pw]]$optimism_mean         <- mean(optimism, na.rm = TRUE)
  elastic_net_results[[pw]]$optimism_n_valid      <- sum(!is.na(optimism))

  cat(sprintf("  %s: R2_apparent=%.3f, optimism=%.3f, R2_corrected=%.3f (%d/%d valid)\n",
              pw, apparent_R2, mean(optimism, na.rm = TRUE), corrected_R2,
              sum(!is.na(optimism)), B_OPT))
}

# ---- Permutation test for Elastic Net ----
cat("\n--- Permutation test for Elastic Net ---\n")
N_PERM_EN <- 1000

for (pw in names(elastic_net_results)) {
  # Per-pathway seed for reproducibility
  pw_idx_en <- which(names(elastic_net_results) == pw)
  set.seed(1000 + pw_idx_en)

  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)

  y <- patient_level[[z_pathway]]
  X <- as.matrix(patient_level[, z_features, drop = FALSE])
  N <- length(y)

  observed_R2 <- elastic_net_results[[pw]]$R2_cv

  # Short-circuit: if observed R2 <= 0, permutation p-value is trivially 1.0
  if (observed_R2 <= 0) {
    cat(sprintf("  %s: R2_cv=%.3f <= 0, skipping permutation (perm_p=1.000)\n",
                pw, observed_R2))
    elastic_net_results[[pw]]$perm_p   <- 1.0
    elastic_net_results[[pw]]$perm_R2s <- rep(NA_real_, N_PERM_EN)
    next
  }

  perm_R2 <- rep(NA_real_, N_PERM_EN)

  for (perm in 1:N_PERM_EN) {
    y_perm <- sample(y)  # shuffle outcome across patients
    perm_preds <- numeric(N)

    for (i in 1:N) {
      # Use same alpha grid search as main LOPO to ensure valid null distribution
      perm_best_alpha  <- 0.5
      perm_best_lambda <- NULL
      perm_best_cvm    <- Inf

      for (alpha_val in seq(0.1, 1.0, by = 0.1)) {
        cv_fit <- tryCatch({
          suppressWarnings(cv.glmnet(X[-i, , drop = FALSE], y_perm[-i],
                    alpha = alpha_val, nfolds = min(5, N - 1)))
        }, error = function(e) NULL)

        if (!is.null(cv_fit)) {
          idx_1se_p <- which(cv_fit$lambda == cv_fit$lambda.1se)
          cvm_1se_p <- cv_fit$cvm[idx_1se_p]
          if (length(cvm_1se_p) > 0 && cvm_1se_p < perm_best_cvm) {
            perm_best_cvm    <- cvm_1se_p
            perm_best_alpha  <- alpha_val
            perm_best_lambda <- cv_fit$lambda.1se
          }
        }
      }

      if (!is.null(perm_best_lambda)) {
        perm_fit <- glmnet(X[-i, , drop = FALSE], y_perm[-i],
                           alpha = perm_best_alpha, lambda = perm_best_lambda)
        perm_preds[i] <- predict(perm_fit, newx = X[i, , drop = FALSE])[1]
      } else {
        perm_preds[i] <- mean(y_perm[-i])
      }
    }

    SS_res_p <- sum((y_perm - perm_preds)^2)
    SS_tot_p <- sum((y_perm - mean(y_perm))^2)
    perm_R2[perm] <- ifelse(SS_tot_p > 0, 1 - SS_res_p / SS_tot_p, NA)
  }

  perm_p <- mean(perm_R2 >= observed_R2, na.rm = TRUE)

  elastic_net_results[[pw]]$perm_p   <- perm_p
  elastic_net_results[[pw]]$perm_R2s <- perm_R2

  cat(sprintf("  %s: R2_cv=%.3f, perm_p=%.4f\n", pw, observed_R2, perm_p))
}

saveRDS(elastic_net_results, file.path(RESULTS_DIR, "elastic_net_results.rds"))

# ============================================================
# STEP 7: SENSITIVITY ANALYSES
# ============================================================

cat("\n========== STEP 7: SENSITIVITY ANALYSES ==========\n")

# --- S1: Alternative zone-to-subcompartment mappings ---
# Implemented as a function to reuse the pipeline

# Load zone-level data for re-aggregation
ssgsea_zone   <- readRDS(file.path(DATA_DIR, "processed", "ssgsea_scores.rds"))
meta_matched  <- readRDS(file.path(DATA_DIR, "processed", "sample_metadata_matched.rds"))
radio_per_sc  <- readRDS(file.path(DATA_DIR, "processed", "radiomics_per_subcompartment.rds"))
pw_lookup     <- readRDS(file.path(DATA_DIR, "processed", "pathway_lookup.rds"))

# Build zone-level ssGSEA data frame (samples x pathways + zone + patient_id)
ssgsea_zone_t <- as.data.frame(t(ssgsea_zone))
ssgsea_zone_t$zone <- meta_matched$zone
ssgsea_zone_t$patient_id <- meta_matched$patient_id

# Helper: re-aggregate ssGSEA with a given zone mapping, merge with radiomics,
# run LMM for specified pathways, return results
run_sensitivity_mapping <- function(zone_map, map_name, ssgsea_zone_df,
                                     radio_per_sc, pw_lookup, lmm_results,
                                     features_pass_cor, agg_fun = mean) {
  cat(sprintf("\n  S1: %s\n", map_name))

  # Aggregate ssGSEA per patient x subcompartment
  agg_list <- list()
  pathway_names_orig <- pw_lookup$pathway_name
  for (sc in names(zone_map)) {
    zones <- zone_map[[sc]]
    sc_data <- ssgsea_zone_df[ssgsea_zone_df$zone %in% zones, ]
    if (nrow(sc_data) == 0) next
    sc_agg <- sc_data %>%
      group_by(patient_id) %>%
      summarise(across(all_of(pathway_names_orig), ~ agg_fun(.x, na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(subcompartment = sc)
    agg_list[[sc]] <- sc_agg
    cat(sprintf("    %s (%s): %d patients\n", sc, paste(zones, collapse="+"), nrow(sc_agg)))
  }
  ssgsea_agg <- bind_rows(agg_list)

  # Rename to pw_N
  ssgsea_for_merge <- ssgsea_agg %>% select(patient_id, subcompartment, all_of(pathway_names_orig))
  colnames(ssgsea_for_merge) <- c("patient_id", "subcompartment", pw_lookup$pw_col)

  # Prepare radiomics in long format (same as script 02)
  matched_pats <- unique(ssgsea_for_merge$patient_id)
  rad_lookup <- readRDS(file.path(DATA_DIR, "processed", "feature_lookup.rds"))
  rad_original <- rad_lookup$feature_name
  rad_new <- rad_lookup$rad_col

  radio_long_list <- list()
  for (sc in names(radio_per_sc)) {
    sc_df <- radio_per_sc[[sc]]
    # Normalize feature names (strip subcompartment tag)
    cn <- colnames(sc_df)
    cn <- sub(paste0("_", sc, "_"), "_", cn)
    colnames(sc_df) <- cn
    sc_df <- sc_df[!duplicated(colnames(sc_df))]
    sc_df$subcompartment <- sc
    sc_df <- sc_df[sc_df$patient_id %in% matched_pats, ]
    radio_long_list[[sc]] <- sc_df
  }
  radio_long <- bind_rows(radio_long_list)
  radio_for_merge <- radio_long %>% select(patient_id, subcompartment, all_of(rad_original))
  colnames(radio_for_merge) <- c("patient_id", "subcompartment", rad_new)

  # Merge
  sens_data <- inner_join(ssgsea_for_merge, radio_for_merge,
                          by = c("patient_id", "subcompartment"))
  sens_data$subcompartment <- factor(sens_data$subcompartment, levels = c("ET", "NET", "ED"))
  cat(sprintf("    Merged: %d rows, %d patients\n", nrow(sens_data),
              length(unique(sens_data$patient_id))))

  if (nrow(sens_data) < 10) {
    cat("    Too few rows for LMM, skipping.\n")
    return(NULL)
  }

  # Z-score within subcompartment
  sens_data <- sens_data %>%
    group_by(subcompartment) %>%
    mutate(across(all_of(intersect(features_pass_cor, colnames(sens_data))),
                  ~ as.numeric(scale(.)), .names = "z_{.col}")) %>%
    ungroup()
  sens_data <- sens_data %>%
    mutate(across(all_of(intersect(pw_lookup$pw_col, colnames(sens_data))),
                  ~ as.numeric(scale(.)), .names = "z_{.col}"))

  # Re-run LMM for pathways that had features in primary analysis
  sens_lmm <- list()
  for (pw in names(lmm_results)) {
    top5 <- lmm_results[[pw]]$feature_names
    z_features <- paste0("z_", top5)
    z_pathway  <- paste0("z_", pw)

    # Check all columns exist
    needed <- c(z_features, z_pathway, "subcompartment", "patient_id")
    if (!all(needed %in% colnames(sens_data))) {
      cat(sprintf("    %s: missing columns, skip\n", pw))
      next
    }

    tryCatch({
      fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
      full_f <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))
      null_f <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

      suppressWarnings(suppressMessages({
        full_ml <- lmer(full_f, data = sens_data, REML = FALSE)
        null_ml <- lmer(null_f, data = sens_data, REML = FALSE)
        lrt <- anova(null_ml, full_ml)

        full_reml <- lmer(full_f, data = sens_data, REML = TRUE)
        r2 <- r2_nakagawa(full_reml, verbose = FALSE)
      }))

      sens_lmm[[pw]] <- list(
        R2_marginal = r2$R2_marginal,
        LRT_p       = lrt$`Pr(>Chisq)`[2]
      )
      cat(sprintf("    %s: R2_m=%.3f, LRT_p=%.4f\n", pw, r2$R2_marginal, lrt$`Pr(>Chisq)`[2]))
    }, error = function(e) {
      cat(sprintf("    %s: FAILED — %s\n", pw, e$message))
    })
  }
  return(sens_lmm)
}

# S1a: ET-only mapping — ET = CT only (drop CTmvp)
s1a_map <- list(ET = c("CT"), NET = c("CTpan"), ED = c("IT", "LE"))
s1a_results <- run_sensitivity_mapping(s1a_map, "S1a: ET=CT only (drop CTmvp)",
                                        ssgsea_zone_t, radio_per_sc, pw_lookup,
                                        lmm_results, features_pass_cor)

# S1b: Conservative mapping — ET = CT, NET = CTpan, ED = LE only (drop IT)
s1b_map <- list(ET = c("CT"), NET = c("CTpan"), ED = c("LE"))
s1b_results <- run_sensitivity_mapping(s1b_map, "S1b: Conservative (ET=CT, ED=LE only)",
                                        ssgsea_zone_t, radio_per_sc, pw_lookup,
                                        lmm_results, features_pass_cor)

saveRDS(list(S1a = s1a_results, S1b = s1b_results),
        file.path(RESULTS_DIR, "sensitivity_s1_alt_mappings.rds"))

# --- S6: Random slopes sensitivity ---
cat("\n--- S6: Random slopes sensitivity ---\n")

for (pw in names(lmm_results)) {
  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)
  fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")

  rs_formula <- as.formula(
    paste(z_pathway, "~", fixed_terms, "+ (1 + subcompartment | patient_id)")
  )

  tryCatch({
    rs_model <- lmer(
      rs_formula,
      data    = analysis_data,
      REML    = TRUE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 20000)
      )
    )

    if (isSingular(rs_model)) {
      cat(sprintf("  %s: Random slopes = SINGULAR (supports simpler model)\n", pw))
    } else {
      rs_r2 <- r2_nakagawa(rs_model, verbose = FALSE)
      ri_r2 <- lmm_results[[pw]]$R2_marginal
      cat(sprintf("  %s: RS R2_m=%.3f vs RI R2_m=%.3f, AIC_RS=%.1f vs AIC_RI=%.1f\n",
                  pw, rs_r2$R2_marginal, ri_r2,
                  AIC(rs_model), lmm_results[[pw]]$AIC))
    }
  }, error = function(e) {
    cat(sprintf("  %s: Random slopes FAILED — %s\n", pw, e$message))
  })
}

# --- S5: Feature selection INSIDE LOPO-CV (bias-free sensitivity) ---
cat("\n--- S5: Feature selection inside LOPO-CV ---\n")

en_results_s5 <- list()

for (pw in target_pathways) {
  z_pathway <- paste0("z_", pw)
  y <- patient_level[[z_pathway]]
  # Use all features that passed correlation filter (not pre-selected top 5)
  all_z_features <- paste0("z_", features_pass_cor)
  X_all <- as.matrix(patient_level[, all_z_features, drop = FALSE])
  N <- length(y)

  predictions_s5 <- numeric(N)

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

    if (length(sig_feats) == 0) {
      predictions_s5[i] <- mean(y_train)
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
      fit_s5 <- glmnet(X_sel_train, y_train, alpha = best_alpha, lambda = best_lambda)
      predictions_s5[i] <- predict(fit_s5, newx = X_sel_test)[1]
    } else {
      predictions_s5[i] <- mean(y_train)
    }
  }

  SS_res_s5 <- sum((y - predictions_s5)^2)
  SS_tot_s5 <- sum((y - mean(y))^2)
  R2_s5 <- 1 - SS_res_s5 / SS_tot_s5

  en_results_s5[[pw]] <- list(R2_cv = R2_s5, predictions = predictions_s5)

  # Compare with pre-selected result
  R2_presel <- elastic_net_results[[pw]]$R2_cv
  cat(sprintf("  %s: R2_inside_cv=%.3f vs R2_preselected=%.3f (delta=%.3f)\n",
              pw, R2_s5, R2_presel, R2_presel - R2_s5))
}

saveRDS(en_results_s5, file.path(RESULTS_DIR, "elastic_net_s5_inside_cv.rds"))

# --- S2: Exclude patients with < 2 subcompartments ---
cat("\n--- S2: Minimum sample count filter ---\n")

# Count subcompartments per patient in the analysis data
sc_per_patient <- analysis_data %>%
  group_by(patient_id) %>%
  summarise(n_sc = n_distinct(subcompartment), .groups = "drop")

cat(sprintf("  Subcompartment distribution: 1 SC=%d, 2 SC=%d, 3 SC=%d patients\n",
            sum(sc_per_patient$n_sc == 1),
            sum(sc_per_patient$n_sc == 2),
            sum(sc_per_patient$n_sc == 3)))

# Exclude patients with only 1 subcompartment (no within-patient variance)
patients_ge2 <- sc_per_patient$patient_id[sc_per_patient$n_sc >= 2]
analysis_data_s2 <- analysis_data[analysis_data$patient_id %in% patients_ge2, ]
cat(sprintf("  After filter (>=2 SC): %d patients, %d rows\n",
            length(patients_ge2), nrow(analysis_data_s2)))

s2_results <- list()
for (pw in names(lmm_results)) {
  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)
  fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
  full_formula <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))
  null_formula <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

  tryCatch({
    suppressWarnings(suppressMessages({
      full_ml <- lmer(full_formula, data = analysis_data_s2, REML = FALSE)
      null_ml <- lmer(null_formula, data = analysis_data_s2, REML = FALSE)
      lrt <- anova(null_ml, full_ml)
      full_reml <- lmer(full_formula, data = analysis_data_s2, REML = TRUE)
      r2 <- r2_nakagawa(full_reml, verbose = FALSE)
    }))

    s2_results[[pw]] <- list(
      R2_marginal = r2$R2_marginal,
      LRT_p       = lrt$`Pr(>Chisq)`[2],
      N_patients  = length(unique(analysis_data_s2$patient_id[
        !is.na(analysis_data_s2[[z_pathway]])]))
    )
    cat(sprintf("  %s: R2_m=%.3f, LRT_p=%.4f (primary: R2_m=%.3f, LRT_p=%.4f)\n",
                pw, r2$R2_marginal, lrt$`Pr(>Chisq)`[2],
                lmm_results[[pw]]$R2_marginal, lmm_results[[pw]]$LRT_p))
  }, error = function(e) {
    cat(sprintf("  %s: FAILED — %s\n", pw, e$message))
  })
}

saveRDS(s2_results, file.path(RESULTS_DIR, "sensitivity_s2_min_samples.rds"))

# --- S3: Median aggregation (instead of mean) ---
cat("\n--- S3: Median aggregation sensitivity ---\n")

# Re-aggregate using median instead of mean via the helper function
s3_map <- list(ET = c("CT", "CTmvp"), NET = c("CTpan"), ED = c("IT", "LE"))
s3_results <- run_sensitivity_mapping(s3_map, "S3: Median aggregation (same mapping)",
                                       ssgsea_zone_t, radio_per_sc, pw_lookup,
                                       lmm_results, features_pass_cor,
                                       agg_fun = median)

saveRDS(s3_results, file.path(RESULTS_DIR, "sensitivity_s3_median_agg.rds"))

# --- S4a-c: Alternative gene set selections ---
cat("\n--- S4a-c: Alternative gene set selections ---\n")

# S4a: Hallmark only (pw_1 through pw_15)
# S4b: Neftel only (pw_16 through pw_19)
# S4c: IvyGAP zone modules only (pw_20 through pw_24)

s4_subsets <- list(
  S4a_Hallmark = paste0("pw_", 1:15),
  S4b_Neftel   = paste0("pw_", 16:19),
  S4c_IvyGAP   = paste0("pw_", 20:24)
)

s4_results <- list()

for (s4_name in names(s4_subsets)) {
  pw_subset <- s4_subsets[[s4_name]]
  # Only re-analyze pathways in this subset that were also in lmm_results
  pw_in_lmm <- intersect(pw_subset, names(lmm_results))
  cat(sprintf("\n  %s: %d pathways in subset, %d with LMM features\n",
              s4_name, length(pw_subset), length(pw_in_lmm)))

  if (length(pw_in_lmm) == 0) {
    cat(sprintf("    No pathways in %s had features — nothing to re-test\n", s4_name))
    s4_results[[s4_name]] <- list(note = "No pathways with features in this subset")
    next
  }

  # Re-apply FDR correction within this subset only
  subset_p <- model_p_values[pw_in_lmm]
  subset_fdr <- p.adjust(subset_p, method = "BH")

  s4_res <- data.frame(
    pathway     = pw_in_lmm,
    R2_marginal = sapply(pw_in_lmm, function(p) lmm_results[[p]]$R2_marginal),
    LRT_p       = subset_p,
    FDR_subset  = subset_fdr,
    FDR_primary = model_fdr[pw_in_lmm],
    stringsAsFactors = FALSE
  )

  n_sig <- sum(subset_fdr < 0.05)
  cat(sprintf("    FDR<0.05 within %s: %d/%d\n", s4_name, n_sig, length(pw_in_lmm)))
  print(s4_res[order(s4_res$FDR_subset), ], row.names = FALSE)

  s4_results[[s4_name]] <- s4_res
}

saveRDS(s4_results, file.path(RESULTS_DIR, "sensitivity_s4_gene_subsets.rds"))

# --- S7: Global standardization (instead of within-subcompartment) ---
cat("\n--- S7: Global standardization sensitivity ---\n")

s7_results <- list()
analysis_data_s7 <- analysis_data  # copy

# Z-score globally instead of within-subcompartment
for (feat in features_pass_cor) {
  vals <- analysis_data_s7[[feat]]
  analysis_data_s7[[paste0("z_", feat)]] <- as.numeric(scale(vals))
}

for (pw in names(lmm_results)) {
  top5 <- lmm_results[[pw]]$feature_names
  z_features <- paste0("z_", top5)
  z_pathway  <- paste0("z_", pw)
  fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
  full_formula <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))
  null_formula <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

  tryCatch({
    full_ml_s7 <- lmer(full_formula, data = analysis_data_s7, REML = FALSE)
    null_ml_s7 <- lmer(null_formula, data = analysis_data_s7, REML = FALSE)
    lrt_s7 <- anova(null_ml_s7, full_ml_s7)

    full_reml_s7 <- lmer(full_formula, data = analysis_data_s7, REML = TRUE)
    r2_s7 <- r2_nakagawa(full_reml_s7, verbose = FALSE)

    s7_results[[pw]] <- list(
      R2_marginal = r2_s7$R2_marginal,
      LRT_p       = lrt_s7$`Pr(>Chisq)`[2]
    )

    cat(sprintf("  %s: S7 R2_m=%.3f (primary=%.3f), S7 LRT_p=%.4f (primary=%.4f)\n",
                pw, r2_s7$R2_marginal, lmm_results[[pw]]$R2_marginal,
                lrt_s7$`Pr(>Chisq)`[2], lmm_results[[pw]]$LRT_p))
  }, error = function(e) {
    cat(sprintf("  %s: S7 FAILED — %s\n", pw, e$message))
  })
}

saveRDS(s7_results, file.path(RESULTS_DIR, "sensitivity_s7_global_std.rds"))

# --- S8 (new): Kenward-Roger df sensitivity ---
cat("\n--- S8: Kenward-Roger df sensitivity ---\n")

s8_results <- list()

if (requireNamespace("pbkrtest", quietly = TRUE)) {
  for (pw in names(lmm_results)) {
    top5 <- lmm_results[[pw]]$feature_names
    z_features <- paste0("z_", top5)
    z_pathway  <- paste0("z_", pw)
    fixed_terms <- paste(c(z_features, "subcompartment"), collapse = " + ")
    full_formula <- as.formula(paste(z_pathway, "~", fixed_terms, "+ (1 | patient_id)"))
    null_formula <- as.formula(paste(z_pathway, "~ subcompartment + (1 | patient_id)"))

    tryCatch({
      full_kr <- lmer(full_formula, data = analysis_data, REML = TRUE)
      null_kr <- lmer(null_formula, data = analysis_data, REML = TRUE)
      kr_test <- pbkrtest::KRmodcomp(full_kr, null_kr)

      s8_results[[pw]] <- list(
        KR_F     = kr_test$stats$Fstat,
        KR_p     = kr_test$stats$p.value,
        LRT_p    = lmm_results[[pw]]$LRT_p
      )

      cat(sprintf("  %s: KR F=%.2f, p=%.4f (vs Satterthwaite LRT p=%.4f)\n",
                  pw, kr_test$stats$Fstat, kr_test$stats$p.value,
                  lmm_results[[pw]]$LRT_p))
    }, error = function(e) {
      cat(sprintf("  %s: KR FAILED — %s\n", pw, e$message))
      s8_results[[pw]] <<- list(KR_F = NA, KR_p = NA, LRT_p = lmm_results[[pw]]$LRT_p)
    })
  }
} else {
  cat("  pbkrtest not installed — skipping Kenward-Roger sensitivity\n")
  cat("  Install with: install.packages('pbkrtest')\n")
}

saveRDS(s8_results, file.path(RESULTS_DIR, "sensitivity_s8_kenward_roger.rds"))

# --- Gene set overlap analysis (Jaccard similarity) ---
cat("\n--- Gene set overlap analysis ---\n")

compute_geneset_overlap <- function(gene_sets) {
  # gene_sets: named list of character vectors (gene symbols)
  n <- length(gene_sets)
  nms <- names(gene_sets)
  jaccard <- matrix(0, n, n, dimnames = list(nms, nms))

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      inter <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
      union <- length(union(gene_sets[[i]], gene_sets[[j]]))
      jaccard[i, j] <- jaccard[j, i] <- ifelse(union > 0, inter / union, 0)
    }
  }
  diag(jaccard) <- 1

  # Identify clusters at threshold 0.30
  dist_mat <- as.dist(1 - jaccard)
  hc <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.70)  # 1 - 0.30 = 0.70 distance threshold

  list(jaccard = jaccard, clusters = clusters, n_independent = length(unique(clusters)))
}

# Execute Jaccard overlap
gene_sets_file <- file.path(DATA_DIR, "processed", "gene_sets_used.rds")
if (file.exists(gene_sets_file)) {
  gene_sets_all <- readRDS(gene_sets_file)
  jaccard_result <- compute_geneset_overlap(gene_sets_all)
  cat(sprintf("  Computed Jaccard overlap for %d gene sets\n", length(gene_sets_all)))
  cat(sprintf("  Independent signal clusters (Jaccard > 0.30): %d\n",
              jaccard_result$n_independent))

  # Save Jaccard matrix as CSV
  jaccard_df <- as.data.frame(round(jaccard_result$jaccard, 4))
  write.csv(jaccard_df, file.path(RESULTS_DIR, "geneset_jaccard_overlap.csv"))
  cat(sprintf("  Jaccard matrix saved to: %s\n",
              file.path(RESULTS_DIR, "geneset_jaccard_overlap.csv")))

  # Report any high-overlap pairs (Jaccard > 0.10)
  jm <- jaccard_result$jaccard
  high_pairs <- which(jm > 0.10 & upper.tri(jm), arr.ind = TRUE)
  if (nrow(high_pairs) > 0) {
    cat("  Gene set pairs with Jaccard > 0.10:\n")
    for (k in seq_len(nrow(high_pairs))) {
      i <- high_pairs[k, 1]; j <- high_pairs[k, 2]
      cat(sprintf("    %s — %s: J=%.3f\n",
                  rownames(jm)[i], colnames(jm)[j], jm[i, j]))
    }
  } else {
    cat("  No gene set pairs with Jaccard > 0.10 (all sets independent)\n")
  }

  saveRDS(jaccard_result, file.path(RESULTS_DIR, "geneset_jaccard_full.rds"))
} else {
  cat("  WARNING: gene_sets_used.rds not found — skipping Jaccard overlap\n")
}

# ============================================================
# STEP 8: TIER CLASSIFICATION
# ============================================================

cat("\n========== STEP 8: RESULT TIER CLASSIFICATION ==========\n")

n_sig_fdr005 <- sum(results_summary$p_FDR < 0.05, na.rm = TRUE)
n_sig_fdr001 <- sum(results_summary$p_FDR < 0.01, na.rm = TRUE)
max_R2m <- max(results_summary$R2_marginal, na.rm = TRUE)
all_fdr_above_010 <- all(results_summary$p_FDR > 0.10, na.rm = TRUE)
all_R2m_below_005 <- all(results_summary$R2_marginal < 0.05, na.rm = TRUE)

# Check permutation results (use LRT chi-squared as primary)
all_perm_above_010 <- all(
  sapply(permutation_results, function(x) x$perm_p_chisq) > 0.10,
  na.rm = TRUE
)
any_perm_below_005 <- any(
  sapply(permutation_results, function(x) x$perm_p_chisq) < 0.05,
  na.rm = TRUE
)

# Check Elastic Net results
best_EN_R2 <- if (length(elastic_net_results) > 0) {
  max(sapply(elastic_net_results, function(x) x$R2_cv), na.rm = TRUE)
} else {
  -1
}

# Classify — Tier 2 requires at least ONE significance indicator (FDR or permutation),
# not just elevated R2m alone (which could be noise without formal significance)
tier <- if (all_fdr_above_010 && all_R2m_below_005 && all_perm_above_010 && best_EN_R2 < 0) {
  "TIER 1: DEFINITIVE NEGATIVE"
} else if (n_sig_fdr005 >= 4 || (max_R2m > 0.15 && any_perm_below_005) || best_EN_R2 > 0.10) {
  "TIER 3: MODERATE POSITIVE"
} else if (n_sig_fdr005 >= 1 || (any_perm_below_005 && max_R2m >= 0.05) || best_EN_R2 >= 0) {
  "TIER 2: WEAK/MARGINAL POSITIVE"
} else {
  "TIER 1: DEFINITIVE NEGATIVE"
}

# Check for Tier 4
if (n_sig_fdr001 >= 8 && max_R2m > 0.25 && best_EN_R2 > 0.20) {
  tier <- "TIER 4: STRONG POSITIVE"
}

cat(sprintf("\n*** RESULT CLASSIFICATION: %s ***\n", tier))
cat(sprintf("  Pathways FDR<0.05: %d\n", n_sig_fdr005))
cat(sprintf("  Pathways FDR<0.01: %d\n", n_sig_fdr001))
cat(sprintf("  Max marginal R2: %.3f\n", max_R2m))
cat(sprintf("  Best Elastic Net R2_cv: %.3f\n", best_EN_R2))

# ============================================================
# SAVE SESSION
# ============================================================

cat("\n========== SAVING SESSION ==========\n")

save.image(file.path(RESULTS_DIR, "03_analysis_session.RData"))
cat(sprintf("Session saved to: %s\n", file.path(RESULTS_DIR, "03_analysis_session.RData")))

# Log session info for reproducibility
sink(file.path(RESULTS_DIR, "03_session_info.txt"))
cat("Session info for 03_statistical_analysis.R\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Seed: 42\n\n"))
sessionInfo()
sink()
cat(sprintf("Session info saved to: %s\n",
            file.path(RESULTS_DIR, "03_session_info.txt")))

cat("\n=== ANALYSIS PIPELINE COMPLETE ===\n")
