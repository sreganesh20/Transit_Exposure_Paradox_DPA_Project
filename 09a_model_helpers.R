# Shared model utilities; sourcing does not read data or train models.

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ranger)
})
SEED <- 5712026
N_FOLDS <- 5L

prepare_features <- function(d) {
  d <- copy(as.data.table(d))
required_cols <- c(
  "query_id",
  "split",
  "spatial_block",
  "coverage_status",
  "eligible_full_coverage",
  "area_km2",
  "center_lat",
  "center_lon",
  "stations_inside_query",
  "total_cta_entries_at_stations_inside_query_2024",
  "crimes_2024"
)
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
stopifnot(
  !anyNA(d$query_id),
  anyDuplicated(d$query_id) == 0,
  !anyNA(d$split),
  !anyNA(d$eligible_full_coverage)
)
numeric_cols <- c(
  "area_km2",
  "center_lat",
  "center_lon",
  "stations_inside_query",
  "total_cta_entries_at_stations_inside_query_2024",
  "crimes_2024"
)
for (col in numeric_cols) {
  if (!is.numeric(d[[col]]) ||
      any(!is.finite(d[[col]]))) {
    stop("Invalid numeric values in: ", col)
  }
}
stopifnot(
  all(d$area_km2 > 0),
  all(d$crimes_2024 >= 0),
  all(d$stations_inside_query >= 0),
  all(d$total_cta_entries_at_stations_inside_query_2024 >= 0)
)
d[, log_area := log(area_km2)]
d[, lat_c := center_lat - 41.88]
d[, lon_c := center_lon + 87.63]
d[, log_stations := log1p(stations_inside_query)]
d[, log_ridership := log1p(
  total_cta_entries_at_stations_inside_query_2024
)]
d[, log_crimes := log1p(crimes_2024)]
  d
}

baseline_features <- c("log_area", "lat_c", "lon_c")
augmented_features <- c(baseline_features, "log_stations", "log_ridership")
target <- "log_crimes"

# Mean Absolute Error
mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}
# Root Mean Squared Error
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}
# Coefficient of determination
r_squared <- function(actual, predicted) {
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  if (sst == 0) {
    return(NA_real_)
  }
  1 - sse / sst
}
# Convert log-scale predictions to crime counts
# Prevent negative predicted crime counts
clip_log_predictions <- function(pred_log) {
  stopifnot(all(is.finite(pred_log)))
  pmax(pred_log, 0)
}

backtransform <- function(pred_log) {
  expm1(clip_log_predictions(pred_log))
}
# Standardized model evaluation
evaluate_predictions <- function(actual_count, pred_log) {
  stopifnot(
    length(actual_count) == length(pred_log),
    length(actual_count) > 0,
    all(is.finite(actual_count)),
    all(is.finite(pred_log)),
    all(actual_count >= 0)
  )
  pred_log <- clip_log_predictions(pred_log)
  actual_log <- log1p(actual_count)
  pred_count <- backtransform(pred_log)
  stopifnot(all(is.finite(pred_count)))
  data.table(
    n = length(actual_count),
    log_mae = mae(actual_log, pred_log),
    log_rmse = rmse(actual_log, pred_log),
    log_r2 = r_squared(actual_log, pred_log),
    count_mae = mae(actual_count, pred_count),
    count_rmse = rmse(actual_count, pred_count),
    count_r2 = r_squared(actual_count, pred_count)
  )
}

run_spatial_cv <- function(
    data,
    model_name,
    features,
    fit_function,
    predict_function
) {
  oof_predictions <- rep(NA_real_, nrow(data))
  prediction_counts <- integer(nrow(data))
  fold_results <- vector("list", N_FOLDS)
  for (f in seq_len(N_FOLDS)) {
    train_idx <- which(data$fold_id != f)
    valid_idx <- which(data$fold_id == f)
    cv_train <- data[train_idx]
    cv_valid <- data[valid_idx]
    # Restrict model inputs to approved predictors
    # and the training target
    fit_data <- as.data.frame(cv_train[, c(features, target), with = FALSE])
    predict_data <- as.data.frame(cv_valid[, features, with = FALSE])
    # Train model using current training folds
    fitted_model <- fit_function(fit_data, features)
    # Predict held-out validation fold
    pred_log <- as.numeric(predict_function(fitted_model, predict_data))
    stopifnot(length(pred_log) == length(valid_idx), all(is.finite(pred_log)))
    pred_log <- clip_log_predictions(pred_log)
    # Store out-of-fold predictions
    oof_predictions[valid_idx] <- pred_log
    prediction_counts[valid_idx] <-
      prediction_counts[valid_idx] + 1L
    # Calculate fold-level metrics
    metrics <- evaluate_predictions(cv_valid$crimes_2024, pred_log)
    metrics[, model := model_name]
    metrics[, fold_id := f]
    fold_results[[f]] <- metrics
  }
  # Every training observation must receive
  # exactly one out-of-fold prediction
  stopifnot(all(prediction_counts == 1L), all(is.finite(oof_predictions)))
  # Pooled cross-validation metrics
  pooled_metrics <- evaluate_predictions(data$crimes_2024, oof_predictions)
  pooled_metrics[, model := model_name]
  # Preserve query-level predictions
  predictions <- data.table(
    query_id = data$query_id,
    spatial_block = data$spatial_block,
    fold_id = data$fold_id,
    model = model_name,
    actual_count = data$crimes_2024,
    actual_log = data$log_crimes,
    predicted_log = oof_predictions,
    predicted_count = backtransform(oof_predictions)
  )
  list(
    pooled_metrics = pooled_metrics,
    fold_metrics = rbindlist(fold_results),
    predictions = predictions
  )
}
fit_lm <- function(data, features) {
  model_formula <- reformulate(termlabels = features, response = target)
  model <- lm(formula = model_formula, data = data)
  if (model$rank < length(coef(model))) {
    stop("Rank-deficient linear regression model.")
  }
  model
}
predict_lm <- function(model, newdata) {
  predict(model, newdata = newdata)
}

GAM_K_GRID <- c(5L, 10L, 15L)
GAM_ENGINE <- Sys.getenv("GAM_ENGINE", unset = "auto")
GAM_BAM_THRESHOLD <- as.integer(Sys.getenv("GAM_BAM_THRESHOLD", unset = "10000"))
GAM_NUM_THREADS <- as.integer(Sys.getenv("GAM_NUM_THREADS", unset = "1"))
stopifnot(
  GAM_ENGINE %in% c("auto", "gam", "bam"),
  !is.na(GAM_BAM_THRESHOLD), GAM_BAM_THRESHOLD >= 1L,
  !is.na(GAM_NUM_THREADS), GAM_NUM_THREADS >= 1L
)
fit_gam <- function(data, features, k) {
  smooth_terms <- vapply(
    features,
    function(feature) {
      n_unique <- length(unique(data[[feature]]))
      # A constant predictor provides no information
      if (n_unique <= 1L) {
        return("")
      }
      # Use a linear term when there are too few
      # distinct values for a smooth
      if (n_unique < 3L) {
        return(feature)
      }
      # Restrict basis dimension to available values
      effective_k <- min(k, n_unique)
      sprintf("s(%s, k = %d, bs = 'tp')", feature, effective_k)
    },
    FUN.VALUE = character(1)
  )
  smooth_terms <- smooth_terms[
    nzchar(smooth_terms)
  ]
  model_formula <- reformulate(termlabels = smooth_terms, response = target)
  use_bam <- GAM_ENGINE == "bam" ||
    (GAM_ENGINE == "auto" && nrow(data) >= GAM_BAM_THRESHOLD)

  if (use_bam) {
    model <- mgcv::bam(
      formula = model_formula,
      data = data,
      family = gaussian(),
      method = "fREML",
      discrete = TRUE,
      nthreads = GAM_NUM_THREADS
    )
  } else {
    model <- mgcv::gam(
      formula = model_formula,
      data = data,
      family = gaussian(),
      method = "REML"
    )
  }
  if (!isTRUE(model$converged)) {
    stop("GAM model did not converge.")
  }
  if (!all(is.finite(coef(model)))) {
    stop("GAM contains non-finite coefficients.")
  }
  if (!all(is.finite(fitted(model)))) {
    stop("GAM contains non-finite fitted values.")
  }
  model
}
predict_gam <- function(model, newdata) {
  predict(model, newdata = newdata, type = "response")
}
tune_gam <- function(data, model_name, features) {
  candidate_results <- vector("list", length(GAM_K_GRID))
  tuning_rows <- vector("list", length(GAM_K_GRID))
  for (i in seq_along(GAM_K_GRID)) {
    current_k <- GAM_K_GRID[i]
    cv_result <- run_spatial_cv(
      data = data,
      model_name = model_name,
      features = features,
      fit_function = function(fit_data, fit_features) {
        fit_gam(data = fit_data, features = fit_features, k = current_k)
      },
      predict_function = predict_gam
    )
    candidate_results[[i]] <- cv_result
    metrics <- copy(cv_result$pooled_metrics)
    metrics[, k := current_k]
    tuning_rows[[i]] <- metrics
  }
  tuning_table <- rbindlist(tuning_rows)
  # Select lowest log MAE.
  # Prefer smaller k in case of an exact tie.
  best_idx <- order(tuning_table$log_mae, tuning_table$k)[1]
  best_k <- tuning_table$k[best_idx]
  best_cv <- candidate_results[[best_idx]]
  setorder(tuning_table, log_mae, k)
  list(best_k = best_k, tuning_table = tuning_table, best_cv = best_cv)
}

RF_NUM_TREES <- 300L
RF_NUM_THREADS <- suppressWarnings(as.integer(
  Sys.getenv("RF_NUM_THREADS", unset = "1")
))
stopifnot(!is.na(RF_NUM_THREADS), RF_NUM_THREADS >= 1L)
RF_NODE_GRID <- c(5L, 15L, 30L)
RF_BASELINE_MTRY <- c(1L, 2L, 3L)
RF_AUGMENTED_MTRY <- c(2L, 3L, 5L)
fit_rf <- function(data, features, mtry, min_node_size) {
  stopifnot(mtry >= 1L, mtry <= length(features), min_node_size >= 1L)
  ranger::ranger(
    formula = reformulate(features, response = target),
    data = data,
    num.trees = RF_NUM_TREES,
    mtry = mtry,
    min.node.size = min_node_size,
    splitrule = "variance",
    replace = TRUE,
    sample.fraction = 1,
    importance = "none",
    write.forest = TRUE,
    num.threads = RF_NUM_THREADS,
    seed = SEED
  )
}
predict_rf <- function(model, newdata) {
  predict(model, data = newdata, num.threads = RF_NUM_THREADS)$predictions
}
tune_rf <- function(data, model_name, features, mtry_grid) {
  tuning_grid <- CJ(mtry = mtry_grid, min_node_size = RF_NODE_GRID)
  candidate_results <- vector("list", nrow(tuning_grid))
  tuning_rows <- vector("list", nrow(tuning_grid))
  for (i in seq_len(nrow(tuning_grid))) {
    current_mtry <- tuning_grid$mtry[i]
    current_node <- tuning_grid$min_node_size[i]
    cv_result <- run_spatial_cv(
      data = data,
      model_name = model_name,
      features = features,
      fit_function = function(fit_data, fit_features) {
        fit_rf(
          data = fit_data,
          features = fit_features,
          mtry = current_mtry,
          min_node_size = current_node
        )
      },
      predict_function = predict_rf
    )
    candidate_results[[i]] <- cv_result
    metrics <- copy(cv_result$pooled_metrics)
    metrics[, `:=`(
      mtry = current_mtry,
      min_node_size = current_node,
      num_trees = RF_NUM_TREES
    )]
    tuning_rows[[i]] <- metrics
  }
  tuning_table <- rbindlist(tuning_rows)
  # Select lowest log MAE.
  # Exact ties favor larger nodes and smaller mtry.
  best_idx <- order(tuning_table$log_mae, -tuning_table$min_node_size, tuning_table$mtry)[1]
  best_mtry <- tuning_table$mtry[best_idx]
  best_node <- tuning_table$min_node_size[best_idx]
  best_cv <- candidate_results[[best_idx]]
  setorder(tuning_table, log_mae, -min_node_size, mtry)
  list(best_mtry = best_mtry, best_node = best_node, tuning_table = tuning_table, best_cv = best_cv)
}
