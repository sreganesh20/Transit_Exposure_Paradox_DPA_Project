# 09a_model_training.R — Transit Exposure Paradox / CSP 571
# Coverage-aware ML framework; log MAE is the primary selection metric.

# 1. Configuration, input and feature preparation ----
cat("\n=== CONFIGURATION, INPUT AND FEATURE PREPARATION ===\n")
source("09a_model_helpers.R")
set.seed(SEED)
INPUT_FILE <- Sys.getenv(
  "ML_INPUT_FILE",
  unset = "data/processed/07a_coverage_aware_workloads_2024.csv"
)
if (!file.exists(INPUT_FILE)) {
  stop("Input dataset not found: ", INPUT_FILE)
}
d <- fread(INPUT_FILE)
d <- prepare_features(d)
# Explicit predictor whitelist prevents target leakage
train <- d[
  split == "UCI_train" &
    eligible_full_coverage == TRUE
]
holdout <- d[
  split == "UCI_spatial_holdout"
]
cta <- d[
  split == "CTA_transfer"
]
cta_primary <- cta[
  eligible_full_coverage == TRUE
]
stopifnot(
  nrow(train) > 0,
  nrow(holdout) > 0,
  nrow(cta_primary) > 0,
  !anyNA(train$spatial_block),
  !anyNA(holdout$spatial_block),
  length(unique(train$spatial_block)) >= N_FOLDS,
  length(intersect(train$spatial_block, holdout$spatial_block)) == 0
)


# 2. Spatial cross-validation folds ----
cat("\n=== SPATIAL CROSS-VALIDATION FOLDS ===\n")
# Count eligible training queries per spatial block
block_counts <- train[
  , .(n_queries = .N),
  by = spatial_block
]
# Deterministic ordering: largest blocks first
setorder(block_counts, -n_queries, spatial_block)
# Initialize fold assignments
block_counts[, fold_id := NA_integer_]
fold_load <- integer(N_FOLDS)
# Assign each complete block to the smallest fold
for (i in seq_len(nrow(block_counts))) {
  selected_fold <- which.min(fold_load)
  block_counts$fold_id[i] <- selected_fold
  fold_load[selected_fold] <-
    fold_load[selected_fold] +
    block_counts$n_queries[i]
}
# Attach fold assignments to training queries only
train[
  block_counts,
  on = "spatial_block",
  fold_id := i.fold_id
]
# Validate fold assignments
stopifnot(!anyNA(train$fold_id), uniqueN(train$fold_id) == N_FOLDS, sum(fold_load) == nrow(train))
# Ensure every spatial block belongs to exactly one fold
block_check <- train[
  , .(n_folds = uniqueN(fold_id)),
  by = spatial_block
]
stopifnot(all(block_check$n_folds == 1L))
# Summarize fold sizes
fold_summary <- train[
  , .(validation_queries = .N, validation_blocks = uniqueN(spatial_block)),
  by = fold_id
]
fold_summary[
  , training_queries :=
    nrow(train) - validation_queries
]
setorder(fold_summary, fold_id)
# Verify geographic separation within every fold
for (f in seq_len(N_FOLDS)) {
  cv_train <- train[fold_id != f]
  cv_valid <- train[fold_id == f]
  shared_blocks <- intersect(cv_train$spatial_block, cv_valid$spatial_block)
  stopifnot(
    nrow(cv_train) > 0,
    nrow(cv_valid) > 0,
    length(shared_blocks) == 0,
    nrow(cv_train) + nrow(cv_valid) == nrow(train)
  )
}
print(fold_summary)

# 3. Evaluation and sanity checks ----
cat("\n=== EVALUATION AND SANITY CHECKS ===\n")
test_actual <- c(0, 10, 100, 1000)
test_pred_log <- log1p(test_actual)
test_metrics <- evaluate_predictions(test_actual, test_pred_log)
stopifnot(
  test_metrics$log_mae < 1e-10,
  test_metrics$log_rmse < 1e-10,
  abs(test_metrics$log_r2 - 1) < 1e-10,
  test_metrics$count_mae < 1e-10,
  test_metrics$count_rmse < 1e-10,
  abs(test_metrics$count_r2 - 1) < 1e-10
)
# Test non-negative count predictions
negative_test <- backtransform(c(-2, -1, 0, 1))
stopifnot(all(negative_test >= 0), all(negative_test[1:3] == 0))

# Verify consistent log/count clipping
clipping_test <- evaluate_predictions(
  actual_count = c(0, 0, 10),
  pred_log = c(-2, -1, log1p(10))
)
stopifnot(
  clipping_test$log_mae < 1e-10,
  clipping_test$log_rmse < 1e-10,
  clipping_test$count_mae < 1e-10,
  clipping_test$count_rmse < 1e-10
)
# Test imperfect predictions
imperfect_test <- evaluate_predictions(actual_count = c(0, 10, 100), pred_log = log1p(c(1, 8, 90)))
stopifnot(
  imperfect_test$log_mae > 0,
  imperfect_test$count_mae > 0,
  imperfect_test$log_r2 < 1,
  imperfect_test$count_r2 < 1
)

# 4. Cross-validation engine and linear models ----
cat("\n=== CROSS-VALIDATION ENGINE AND LINEAR MODELS ===\n")
lm_baseline_cv <- run_spatial_cv(
  data = train,
  model_name = "LM-B",
  features = baseline_features,
  fit_function = fit_lm,
  predict_function = predict_lm
)
lm_augmented_cv <- run_spatial_cv(
  data = train,
  model_name = "LM-A",
  features = augmented_features,
  fit_function = fit_lm,
  predict_function = predict_lm
)
lm_log_mae_improvement <- 100 * (
  lm_baseline_cv$pooled_metrics$log_mae -
    lm_augmented_cv$pooled_metrics$log_mae
) / lm_baseline_cv$pooled_metrics$log_mae
# Validate out-of-fold prediction coverage
stopifnot(
  nrow(lm_baseline_cv$predictions) == nrow(train),
  nrow(lm_augmented_cv$predictions) == nrow(train),
  identical(lm_baseline_cv$predictions$query_id, lm_augmented_cv$predictions$query_id),
  identical(lm_baseline_cv$predictions$fold_id, lm_augmented_cv$predictions$fold_id)
)

# 5. GAM training and tuning ----
cat("\n=== GAM TRAINING AND TUNING ===\n")
gam_baseline <- tune_gam(data = train, model_name = "GAM-B", features = baseline_features)
gam_baseline_cv <- gam_baseline$best_cv
gam_augmented <- tune_gam(data = train, model_name = "GAM-A", features = augmented_features)
gam_augmented_cv <- gam_augmented$best_cv
gam_log_mae_improvement <- 100 * (
  gam_baseline_cv$pooled_metrics$log_mae -
    gam_augmented_cv$pooled_metrics$log_mae
) / gam_baseline_cv$pooled_metrics$log_mae
stopifnot(
  gam_baseline$best_k %in% GAM_K_GRID,
  gam_augmented$best_k %in% GAM_K_GRID,
  nrow(gam_baseline$tuning_table) ==
    length(GAM_K_GRID),
  nrow(gam_augmented$tuning_table) ==
    length(GAM_K_GRID),
  nrow(gam_baseline_cv$predictions) ==
    nrow(train),
  nrow(gam_augmented_cv$predictions) ==
    nrow(train),
  identical(gam_baseline_cv$predictions$query_id, lm_baseline_cv$predictions$query_id),
  identical(gam_augmented_cv$predictions$fold_id, lm_augmented_cv$predictions$fold_id)
)


# 6. Random forest training and tuning ----
cat("\n=== RANDOM FOREST TRAINING AND TUNING ===\n")
rf_baseline <- tune_rf(
  data = train,
  model_name = "RF-B",
  features = baseline_features,
  mtry_grid = RF_BASELINE_MTRY
)
rf_baseline_cv <- rf_baseline$best_cv
rf_augmented <- tune_rf(
  data = train,
  model_name = "RF-A",
  features = augmented_features,
  mtry_grid = RF_AUGMENTED_MTRY
)
rf_augmented_cv <- rf_augmented$best_cv
rf_log_mae_improvement <- 100 * (
  rf_baseline_cv$pooled_metrics$log_mae -
    rf_augmented_cv$pooled_metrics$log_mae
) / rf_baseline_cv$pooled_metrics$log_mae
stopifnot(
  rf_baseline$best_mtry %in% RF_BASELINE_MTRY,
  rf_augmented$best_mtry %in% RF_AUGMENTED_MTRY,
  rf_baseline$best_node %in% RF_NODE_GRID,
  rf_augmented$best_node %in% RF_NODE_GRID,
  nrow(rf_baseline$tuning_table) ==
    length(RF_BASELINE_MTRY) * length(RF_NODE_GRID),
  nrow(rf_augmented$tuning_table) ==
    length(RF_AUGMENTED_MTRY) * length(RF_NODE_GRID)
)

# 7. Consolidated cross-validation results ----
cat("\n=== CONSOLIDATED CROSS-VALIDATION RESULTS ===\n")
all_cv_results <- list(
  "LM-B" = lm_baseline_cv,
  "LM-A" = lm_augmented_cv,
  "GAM-B" = gam_baseline_cv,
  "GAM-A" = gam_augmented_cv,
  "RF-B" = rf_baseline_cv,
  "RF-A" = rf_augmented_cv
)
# Validate identical evaluation workloads
for (model_name in names(all_cv_results)) {
  model_predictions <-
    all_cv_results[[model_name]]$predictions
  stopifnot(
    identical(model_predictions$query_id, train$query_id),
    identical(model_predictions$fold_id, train$fold_id),
    nrow(model_predictions) == nrow(train),
    all(is.finite(model_predictions$predicted_log))
  )
}
# Combine existing results without retraining models
cv_comparison <- rbindlist(lapply(all_cv_results, function(x) x$pooled_metrics))
all_fold_metrics <- rbindlist(lapply(all_cv_results, function(x) x$fold_metrics))
all_oof_predictions <- rbindlist(lapply(all_cv_results, function(x) x$predictions))
setcolorder(
  cv_comparison,
  c("model", "n", "log_mae", "log_rmse", "log_r2", "count_mae", "count_rmse", "count_r2")
)
setorder(cv_comparison, log_mae)
# Compare mobility augmentation
mobility_comparison <- data.table(
  model_family = c("Linear Regression", "GAM", "Random Forest"),
  baseline_log_mae = c(
    lm_baseline_cv$pooled_metrics$log_mae,
    gam_baseline_cv$pooled_metrics$log_mae,
    rf_baseline_cv$pooled_metrics$log_mae
  ),
  augmented_log_mae = c(
    lm_augmented_cv$pooled_metrics$log_mae,
    gam_augmented_cv$pooled_metrics$log_mae,
    rf_augmented_cv$pooled_metrics$log_mae
  ),
  improvement_pct = c(lm_log_mae_improvement, gam_log_mae_improvement, rf_log_mae_improvement)
)
# Final cross-validation assertions
stopifnot(
  nrow(cv_comparison) == 6L,
  nrow(all_fold_metrics) == 6L * N_FOLDS,
  nrow(all_oof_predictions) == 6L * nrow(train),
  anyDuplicated(all_oof_predictions[, .(model, query_id)]) == 0L
)
print(cv_comparison)

# 8. Final refit and predictions ----
cat("\n=== FINAL REFIT AND PREDICTIONS ===\n")
baseline_fit_data <- as.data.frame(train[, c(baseline_features, target), with = FALSE])
augmented_fit_data <- as.data.frame(train[, c(augmented_features, target), with = FALSE])
final_models <- list(
  "LM-B" = fit_lm(baseline_fit_data, baseline_features),
  "LM-A" = fit_lm(augmented_fit_data, augmented_features),
  "GAM-B" = fit_gam(baseline_fit_data, baseline_features, k = gam_baseline$best_k),
  "GAM-A" = fit_gam(augmented_fit_data, augmented_features, k = gam_augmented$best_k),
  "RF-B" = fit_rf(
    baseline_fit_data, baseline_features,
    mtry = rf_baseline$best_mtry,
    min_node_size = rf_baseline$best_node
  ),
  "RF-A" = fit_rf(
    augmented_fit_data, augmented_features,
    mtry = rf_augmented$best_mtry,
    min_node_size = rf_augmented$best_node
  )
)
final_fit_summary <- data.table(
  model = names(final_models),
  training_queries = vapply(
    final_models,
    function(m) {
      if (inherits(m, "ranger")) {
        as.integer(m$num.samples)
      } else {
        as.integer(nobs(m))
      }
    },
    integer(1)
  ),
  model_class = vapply(final_models, function(m) class(m)[1], character(1))
)
stopifnot(
  length(final_models) == 6L,
  all(final_fit_summary$training_queries == nrow(train)),
  identical(names(final_models), names(all_cv_results))
)
prediction_functions <- list(
  "LM-B" = predict_lm,
  "LM-A" = predict_lm,
  "GAM-B" = predict_gam,
  "GAM-A" = predict_gam,
  "RF-B" = predict_rf,
  "RF-A" = predict_rf
)
generate_predictions <- function(workload, workload_name) {
  results <- lapply(names(final_models), function(model_name) {
    features <- if (endsWith(model_name, "-B")) {
      baseline_features
    } else {
      augmented_features
    }
    newdata <- as.data.frame(workload[, features, with = FALSE])
    pred_log <- as.numeric(prediction_functions[[model_name]](final_models[[model_name]], newdata))
    stopifnot(length(pred_log) == nrow(workload), all(is.finite(pred_log)))
    pred_log <- clip_log_predictions(pred_log)
    data.table(
      query_id = workload$query_id,
      spatial_block = workload$spatial_block,
      workload = workload_name,
      coverage_status = workload$coverage_status,
      eligible_full_coverage = workload$eligible_full_coverage,
      model = model_name,
      actual_count = workload$crimes_2024,
      actual_log = workload$log_crimes,
      predicted_log = pred_log,
      predicted_count = backtransform(pred_log)
    )
  })
  rbindlist(results)
}
holdout_predictions <- generate_predictions(holdout, "UCI_spatial_holdout")
cta_predictions <- generate_predictions(cta, "CTA_transfer")
final_predictions <- rbindlist(list(holdout_predictions, cta_predictions))
stopifnot(
  nrow(holdout_predictions) == 6L * nrow(holdout),
  nrow(cta_predictions) == 6L * nrow(cta),
  all(is.finite(final_predictions$predicted_count)),
  all(final_predictions$predicted_count >= 0),
  anyDuplicated(final_predictions[, .(model, query_id)]) == 0L
)
prediction_summary <- final_predictions[
  , .(
    total_predictions = .N,
    eligible_predictions = sum(eligible_full_coverage),
    nonnegative_predictions = sum(predicted_count >= 0)
  ),
  by = .(workload, model)
]

# 9. Holdout and transfer evaluation ----
cat("\n=== HOLDOUT AND TRANSFER EVALUATION ===\n")
evaluation_data <- final_predictions[
  eligible_full_coverage == TRUE
]
final_metrics <- evaluation_data[
  , evaluate_predictions(actual_count, predicted_log),
  by = .(workload, model)
]
setorder(final_metrics, workload, log_mae)
stopifnot(
  nrow(final_metrics) == 12L,
  all(final_metrics[ workload == "UCI_spatial_holdout", n ] == sum(holdout$eligible_full_coverage)),
  all(final_metrics[ workload == "CTA_transfer", n ] == nrow(cta_primary))
)

# Previously unused eligible UCI holdout evaluation
unused_holdout_metrics <- NULL
unused_holdout_mobility_comparison <- NULL

if ("is_poc_query" %in% names(holdout)) {
  stopifnot(
    !anyNA(holdout$is_poc_query),
    all(holdout$is_poc_query %in% c(TRUE, FALSE))
  )

  unused_ids <- holdout[
    eligible_full_coverage == TRUE & is_poc_query == FALSE,
    query_id
  ]

  if (length(unused_ids) > 0L) {
    unused_eval <- evaluation_data[
      workload == "UCI_spatial_holdout" &
        query_id %in% unused_ids
    ]

    stopifnot(
      nrow(unused_eval) == length(final_models) * length(unused_ids),
      anyDuplicated(unused_eval[, .(model, query_id)]) == 0L
    )

    unused_holdout_metrics <- unused_eval[
      , evaluate_predictions(actual_count, predicted_log),
      by = .(workload, model)
    ]

    stopifnot(
      nrow(unused_holdout_metrics) == length(final_models),
      all(unused_holdout_metrics$n == length(unused_ids))
    )

    setorder(unused_holdout_metrics, log_mae)

    unused_paired <- copy(unused_holdout_metrics)
    unused_paired[, model_family := sub("-[BA]$", "", model)]
    unused_paired[, configuration := ifelse(
      endsWith(model, "-B"), "baseline", "augmented"
    )]

    unused_holdout_mobility_comparison <- dcast(
      unused_paired,
      workload + model_family ~ configuration,
      value.var = "log_mae"
    )

    unused_holdout_mobility_comparison[, improvement_pct :=
      ifelse(baseline > 0,
        100 * (baseline - augmented) / baseline,
        NA_real_)
    ]

    stopifnot(
      nrow(unused_holdout_mobility_comparison) == 3L
    )

    cat("\nPreviously unused eligible holdout queries:",
        length(unused_ids), "\n")
    print(unused_holdout_metrics)
    print(unused_holdout_mobility_comparison)
  } else {
    cat("\nNo previously unused eligible holdout queries.\n")
  }
} else {
  cat("\nis_poc_query unavailable; unused holdout evaluation skipped.\n")
}

# Compare baseline and augmented performance
paired_metrics <- copy(final_metrics)
paired_metrics[, model_family := sub("-[BA]$", "", model)]
paired_metrics[, configuration := ifelse(
  endsWith(model, "-B"), "baseline", "augmented"
)]
final_mobility_comparison <- dcast(
  paired_metrics,
  workload + model_family ~ configuration,
  value.var = "log_mae"
)
final_mobility_comparison[, improvement_pct :=
                            100 * (baseline - augmented) / baseline
]
stopifnot(
  nrow(final_mobility_comparison) == 6L,
  all(is.finite(final_mobility_comparison$improvement_pct))
)
print(final_metrics)
print(final_mobility_comparison)

# 10. Marginal training-support diagnostics ----
cat("\n=== MARGINAL TRAINING-SUPPORT DIAGNOSTICS ===\n")
support_bounds <- rbindlist(lapply(
  augmented_features,
  function(feature) {
    data.table(
      feature = feature,
      lower = as.numeric(quantile(train[[feature]], 0.01)),
      upper = as.numeric(quantile(train[[feature]], 0.99))
    )
  }
))
diagnose_support <- function(workload, workload_name) {
  result <- data.table(
    query_id = workload$query_id,
    workload = workload_name,
    eligible_full_coverage = workload$eligible_full_coverage
  )
  for (feature in augmented_features) {
    bound_idx <- match(feature, support_bounds$feature)
    stopifnot(length(bound_idx) == 1L, !is.na(bound_idx))
    bounds <- support_bounds[bound_idx, ]
    stopifnot(nrow(bounds) == 1L)
    outside_flag <- workload[[feature]] < bounds$lower |
      workload[[feature]] > bounds$upper
    result[, (paste0(feature, "_outside")) := outside_flag]
  }
  baseline_cols <- paste0(baseline_features, "_outside")
  augmented_cols <- paste0(augmented_features, "_outside")
  result[, baseline_outside :=
           rowSums(.SD) > 0, .SDcols = baseline_cols
  ]
  result[, augmented_outside :=
           rowSums(.SD) > 0, .SDcols = augmented_cols
  ]
  result
}
support_diagnostics <- rbindlist(list(
  diagnose_support(holdout, "UCI_spatial_holdout"),
  diagnose_support(cta, "CTA_transfer")
))
support_summary <- support_diagnostics[
  eligible_full_coverage == TRUE,
  .(
    n_queries = .N,
    baseline_outside_n = sum(baseline_outside),
    baseline_outside_pct = 100 * mean(baseline_outside),
    augmented_outside_n = sum(augmented_outside),
    augmented_outside_pct = 100 * mean(augmented_outside)
  ),
  by = workload
]
stopifnot(
  nrow(support_diagnostics) == nrow(holdout) + nrow(cta),
  !anyNA(support_diagnostics),
  all(!support_diagnostics$baseline_outside | support_diagnostics$augmented_outside)
)
support_eval <- merge(
  evaluation_data,
  support_diagnostics[, .(query_id, augmented_outside)],
  by = "query_id",
  all.x = TRUE
)
stopifnot(nrow(support_eval) == nrow(evaluation_data), !anyNA(support_eval$augmented_outside))
support_eval[, support_group := ifelse(
  augmented_outside,
  "Outside bounds",
  "Inside bounds"
)]
support_metrics <- support_eval[
  , evaluate_predictions(actual_count, predicted_log),
  by = .(workload, support_group, model)
]
setorder(support_metrics, workload, support_group, log_mae)
stopifnot(nrow(support_metrics) == uniqueN(support_eval[, .(workload, support_group, model)]), sum(support_metrics$n) == nrow(evaluation_data))
print(support_summary)

# 11. Paired geographic block evaluation and exploratory bootstrap ----
cat("\n=== PAIRED GEOGRAPHIC BLOCK EVALUATION AND EXPLORATORY BOOTSTRAP ===\n")
# Reuse the original Step 5 geographic grid
tile_id <- function(x) {
  paste(floor((x$center_lat - 41.0) / 0.020), floor((x$center_lon + 88.0) / 0.025), sep = "_")
}
# Verify reproduction of original UCI blocks
stopifnot(
  all(tile_id(train) == train$spatial_block),
  all(tile_id(holdout) == holdout$spatial_block)
)
# Assign CTA blocks for evaluation only
cta_block_map <- setNames(tile_id(cta), cta$query_id)
paired_errors <- copy(evaluation_data)
paired_errors[workload == "CTA_transfer",
              spatial_block := unname(cta_block_map[query_id])
]
stopifnot(!anyNA(paired_errors$spatial_block))
paired_errors[, absolute_log_error :=
                abs(actual_log - predicted_log)
]
paired_errors[, model_family :=
                sub("-[BA]$", "", model)
]
paired_errors[, configuration := ifelse(
  endsWith(model, "-B"), "baseline", "augmented"
)]
paired_errors <- dcast(
  paired_errors,
  workload + spatial_block + query_id + model_family ~
    configuration,
  value.var = "absolute_log_error"
)
stopifnot(
  nrow(paired_errors) == 3L * (sum(holdout$eligible_full_coverage) + nrow(cta_primary)),
  all(is.finite(paired_errors$baseline)),
  all(is.finite(paired_errors$augmented))
)
paired_errors[, error_improvement :=
                baseline - augmented
]
block_comparison <- paired_errors[
  , .(
    n_queries = .N,
    baseline_mae = mean(baseline),
    augmented_mae = mean(augmented),
    mean_improvement = mean(error_improvement)
  ),
  by = .(workload, model_family, spatial_block)
]
block_summary <- block_comparison[
  , .(
    n_blocks = .N,
    n_queries = sum(n_queries),
    improved_blocks = sum(mean_improvement > 0),
    worsened_blocks = sum(mean_improvement < 0),
    tied_blocks = sum(mean_improvement == 0),
    mean_block_improvement = mean(mean_improvement),
    query_weighted_improvement = weighted.mean(mean_improvement, n_queries)
  ),
  by = .(workload, model_family)
]
stopifnot(
  nrow(block_summary) == 6L,
  all(
    block_summary$improved_blocks +
      block_summary$worsened_blocks +
      block_summary$tied_blocks ==
      block_summary$n_blocks
  )
)
BOOT_REPS <- 5000L
set.seed(SEED)
bootstrap_results <- block_comparison[
  , {
    block_n <- n_queries
    block_gain <- mean_improvement
    k <- .N
    boot_gain <- replicate(BOOT_REPS, {
      idx <- sample.int(k, k, replace = TRUE)
      weighted.mean(block_gain[idx], block_n[idx])
    })
    .(
      n_blocks = k,
      observed_improvement = weighted.mean(block_gain, block_n),
      ci_lower = as.numeric(quantile(boot_gain, 0.025)),
      ci_upper = as.numeric(quantile(boot_gain, 0.975))
    )
  },
  by = .(workload, model_family)
]
stopifnot(
  nrow(bootstrap_results) == 6L,
  all(is.finite(bootstrap_results$ci_lower)),
  all(is.finite(bootstrap_results$ci_upper)),
  all(bootstrap_results$ci_lower <= bootstrap_results$ci_upper)
)
print(block_summary)
print(bootstrap_results)

# 12. Export reproducible results ----
cat("\n=== EXPORT REPRODUCIBLE RESULTS ===\n")
OUTPUT_DIR <- Sys.getenv(
  "ML_OUTPUT_DIR",
  unset = "outputs/09a_ml_framework_validation"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
export_tables <- list(
  spatial_folds = block_counts,
  cv_metrics = cv_comparison,
  cv_fold_metrics = all_fold_metrics,
  cv_oof_predictions = all_oof_predictions,
  cv_mobility_comparison = mobility_comparison,
  gam_baseline_tuning = gam_baseline$tuning_table,
  gam_augmented_tuning = gam_augmented$tuning_table,
  rf_baseline_tuning = rf_baseline$tuning_table,
  rf_augmented_tuning = rf_augmented$tuning_table,
  final_fit_summary = final_fit_summary,
  final_predictions = final_predictions,
  final_metrics = final_metrics,
  final_mobility_comparison = final_mobility_comparison,
  support_bounds = support_bounds,
  support_diagnostics = support_diagnostics,
  support_summary = support_summary,
  support_metrics = support_metrics,
  paired_errors = paired_errors,
  block_comparison = block_comparison,
  block_summary = block_summary,
  bootstrap_results = bootstrap_results
)

if (!is.null(unused_holdout_metrics)) {
  export_tables$unused_holdout_ids <- data.table(query_id = unused_ids)
  export_tables$unused_holdout_metrics <-
    unused_holdout_metrics
  export_tables$unused_holdout_mobility_comparison <-
    unused_holdout_mobility_comparison
}

export_tables$selected_hyperparameters <- data.table(
  model = names(final_models),
  k = c(NA_integer_, NA_integer_, gam_baseline$best_k, gam_augmented$best_k, NA_integer_, NA_integer_),
  mtry = c(rep(NA_integer_, 4), rf_baseline$best_mtry, rf_augmented$best_mtry),
  min_node_size = c(rep(NA_integer_, 4), rf_baseline$best_node, rf_augmented$best_node),
  num_trees = c(rep(NA_integer_, 4), RF_NUM_TREES, RF_NUM_TREES),
  gam_engine = c(NA_character_, NA_character_, if (inherits(final_models[["GAM-B"]], "bam")) "bam" else "gam", if (inherits(final_models[["GAM-A"]], "bam")) "bam" else "gam", NA_character_, NA_character_)
)
export_tables$run_metadata <- data.table(
  input_file = INPUT_FILE,
  input_md5 = unname(tools::md5sum(INPUT_FILE)),
  seed = SEED, n_folds = N_FOLDS,
  training_queries = nrow(train),
  gam_engine = GAM_ENGINE, gam_bam_threshold = GAM_BAM_THRESHOLD,
  gam_num_threads = GAM_NUM_THREADS, rf_num_threads = RF_NUM_THREADS
)

export_manifest <- rbindlist(lapply(
  names(export_tables),
  function(table_name) {
    output_file <- file.path(OUTPUT_DIR, paste0(table_name, ".csv"))
    fwrite(export_tables[[table_name]], output_file)
    data.table(
      file = basename(output_file),
      rows = nrow(export_tables[[table_name]]),
      columns = ncol(export_tables[[table_name]]),
      exists = file.exists(output_file)
    )
  }
))
stopifnot(nrow(export_manifest) == length(export_tables), all(export_manifest$exists), all(export_manifest$rows > 0))
fwrite(export_manifest, file.path(OUTPUT_DIR, "export_manifest.csv"))
print(export_manifest)
cat("\nAll checks passed. CSV files exported:", nrow(export_manifest) + 1L, "\n")
