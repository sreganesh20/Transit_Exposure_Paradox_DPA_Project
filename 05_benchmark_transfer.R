# CSP 571 - The Transit Exposure Paradox
# STEP 5: Does a model trained on original UCI spatial queries generalize
# to new CTA-centered rectangular query workloads?
#
# Runs from the project root. Input created by 04b_augment_uci_queries.R.
# Uses base R only (no additional packages, downloads, or paid APIs).
#
# IMPORTANT DESIGN:
#   Outcome: re-computed 2024 crime count for EACH query (same data snapshot).
#   NOT outcome: UCI's original historical 'count'.
#   Training data: UCI rectangles ONLY, with geographically held-out UCI tests.
#   Transfer test: CTA-centered squares, which are NEVER used in training.
#   Prediction models estimate query counts, NOT future crime or passenger risk.

input <- file.path('data', 'processed', '04b_spatial_query_workloads_2024.csv')
output_dir <- 'outputs'
processed_dir <- file.path('data', 'processed')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input)) {
  stop('Missing ', input, '\nRun Step 4B first, from the project root: ', getwd())
}

cat('Loading enriched UCI + CTA rectangular workload...\n')
d <- read.csv(input, stringsAsFactors = FALSE)
needed <- c('query_id', 'query_source', 'area_km2', 'center_lon', 'center_lat',
            'crimes_2024', 'stations_inside_query',
            'total_cta_entries_at_stations_inside_query_2024',
            'station_id', 'query_halfwidth_m')
missing <- setdiff(needed, names(d))
if (length(missing)) stop('Missing columns: ', paste(missing, collapse = ', '))
for (nm in c('area_km2', 'center_lon', 'center_lat', 'crimes_2024',
             'stations_inside_query',
             'total_cta_entries_at_stations_inside_query_2024',
             'query_halfwidth_m')) {
  d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
}
if (anyDuplicated(d$query_id)) stop('Duplicate query_id in Step 4B dataset.')
if (any(!is.finite(d$area_km2)) || any(d$area_km2 <= 0) ||
    any(!is.finite(d$center_lon)) || any(!is.finite(d$center_lat)) ||
    any(!is.finite(d$crimes_2024)) || any(d$crimes_2024 < 0) ||
    any(!is.finite(d$stations_inside_query)) ||
    any(!is.finite(d$total_cta_entries_at_stations_inside_query_2024)) ||
    any(d$stations_inside_query < 0) ||
    any(d$total_cta_entries_at_stations_inside_query_2024 < 0)) {
  stop('Required modeling features contain invalid values.')
}
uci <- d[d$query_source == 'UCI_original_rectangle', , drop = FALSE]
cta <- d[d$query_source == 'CTA_centered_rectangle', , drop = FALSE]
if (nrow(uci) < 250L || nrow(cta) < 100L) {
  stop('Unexpected sizes: UCI=', nrow(uci), ', CTA=', nrow(cta), '. Verify Step 4B.')
}
if (anyNA(cta$station_id) || anyNA(cta$query_halfwidth_m)) {
  stop('Missing CTA station IDs or query halfwidths.')
}
if (!setequal(unique(cta$query_halfwidth_m), c(250, 500, 750))) {
  stop('Expected CTA squares with halfwidths 250, 500 and 750 m.')
}
cat('UCI rectangles: ', nrow(uci), '; CTA squares: ', nrow(cta), '\n', sep = '')

# --- 1. Reproducible spatial-block holdout -----------------------------
# Splitting randomly by query ROW would create heavy leakage because nearby
# UCI rectangles overlap. We instead hold out entire query-center tiles.
# Blocks are approximate geographical cells, NOT fully non-overlapping polygons.
tile_id <- function(x) {
  paste(floor((x$center_lat - 41.0) / 0.020),
        floor((x$center_lon + 88.0) / 0.025), sep = '_')
}
uci$spatial_block <- tile_id(uci)
blocks <- unique(uci$spatial_block)
if (length(blocks) < 8L) {
  stop('Too few geographic blocks for spatial holdout. Blocks: ', length(blocks))
}
set.seed(5712024)
selected_test_blocks <- NULL
# Find a grouping that yields a meaningful test fraction without splitting a tile.
for (attempt in seq_len(250)) {
  proposed <- sample(blocks, max(2L, round(length(blocks) * 0.25)))
  share <- mean(uci$spatial_block %in% proposed)
  if (share >= 0.15 && share <= 0.35) {
    selected_test_blocks <- proposed
    break
  }
}
if (is.null(selected_test_blocks)) {
  stop('Could not create a balanced spatial-block holdout. Inspect spatial blocks.')
}
uci$split <- ifelse(uci$spatial_block %in% selected_test_blocks,
                    'UCI_spatial_holdout', 'UCI_train')
train <- uci[uci$split == 'UCI_train', , drop = FALSE]
holdout <- uci[uci$split == 'UCI_spatial_holdout', , drop = FALSE]
cat('Training UCI queries: ', nrow(train),
    '; geographically held-out UCI queries: ', nrow(holdout), '\n', sep = '')
write.csv(uci[, c('query_id', 'spatial_block', 'split', 'area_km2',
                  'center_lon', 'center_lat')],
          file.path(output_dir, '05_uci_spatial_split.csv'), row.names = FALSE)

# --- 2. Two simple, interpretable models -------------------------------
# Transform the target to reduce count skew, and use ONLY query information
# available before computing the exact query answer.
features <- function(x) {
  x$log_area <- log(x$area_km2)
  x$lat_centered <- x$center_lat - 41.88
  x$lon_centered <- x$center_lon + 87.63
  x$log_station_count <- log1p(x$stations_inside_query)
  x$log_entries <- log1p(x$total_cta_entries_at_stations_inside_query_2024)
  x$log_answer <- log1p(x$crimes_2024)
  x
}
train <- features(train)
holdout <- features(holdout)
cta <- features(cta)

baseline_formula <- log_answer ~ log_area + lat_centered + lon_centered
augmented_formula <- log_answer ~ log_area + lat_centered + lon_centered +
  log_station_count + log_entries
baseline <- lm(baseline_formula, data = train)
augmented <- lm(augmented_formula, data = train)
if (anyNA(coef(baseline)) || anyNA(coef(augmented))) {
  stop('Model coefficients are undefined. Inspect feature collinearity.')
}

# Back-transform predicted log(1 + count). Clamp to nonnegative reasonable
# bounds to avoid negative or numerically overflowing predicted counts.
attach_predictions <- function(x) {
  x$pred_baseline_log <- pmax(0, pmin(log1p(1e8),
                                      as.numeric(predict(baseline, newdata = x))))
  x$pred_augmented_log <- pmax(0, pmin(log1p(1e8),
                                       as.numeric(predict(augmented, newdata = x))))
  if (any(!is.finite(x$pred_baseline_log)) ||
      any(!is.finite(x$pred_augmented_log))) {
    stop('Non-finite predictions encountered.')
  }
  x$pred_baseline_count <- expm1(x$pred_baseline_log)
  x$pred_augmented_count <- expm1(x$pred_augmented_log)
  x$abs_log_error_baseline <- abs(x$log_answer - x$pred_baseline_log)
  x$abs_log_error_augmented <- abs(x$log_answer - x$pred_augmented_log)
  x
}
holdout <- attach_predictions(holdout)
cta <- attach_predictions(cta)

# --- 3. Accuracy by test set -------------------------------------------
score <- function(x, label, method) {
  pred <- if (method == 'Baseline') x$pred_baseline_count else x$pred_augmented_count
  pred_log <- if (method == 'Baseline') x$pred_baseline_log else x$pred_augmented_log
  data.frame(
    evaluation_set = label, model = method, n_queries = nrow(x),
    mae_count = mean(abs(x$crimes_2024 - pred)),
    median_absolute_count_error = median(abs(x$crimes_2024 - pred)),
    rmse_count = sqrt(mean((x$crimes_2024 - pred)^2)),
    mae_log1p_count = mean(abs(log1p(x$crimes_2024) - pred_log)),
    rmsle = sqrt(mean((log1p(x$crimes_2024) - pred_log)^2)),
    stringsAsFactors = FALSE
  )
}
sets <- list(UCI_spatial_holdout = holdout,
             CTA_all = cta,
             CTA_250m_square = cta[cta$query_halfwidth_m == 250, ],
             CTA_500m_square = cta[cta$query_halfwidth_m == 500, ],
             CTA_750m_square = cta[cta$query_halfwidth_m == 750, ])
evaluation <- do.call(rbind, lapply(names(sets), function(name) {
  x <- sets[[name]]
  rbind(score(x, name, 'Baseline'), score(x, name, 'Augmented'))
}))
rownames(evaluation) <- NULL
write.csv(evaluation, file.path(output_dir, '05_evaluation_metrics.csv'),
          row.names = FALSE)

# Area matching: NOTE this is a descriptive diagnostic, not a randomized
# matched experiment. Geography and query aspect ratios still differ.
area_checks <- do.call(rbind, lapply(c(250, 500, 750), function(radius) {
  target_area <- (2 * radius)^2 / 1e6
  matched <- holdout[holdout$area_km2 >= target_area / 1.5 &
                      holdout$area_km2 <= target_area * 1.5, , drop = FALSE]
  target <- cta[cta$query_halfwidth_m == radius, , drop = FALSE]
  data.frame(
    cta_square_halfwidth_m = radius, cta_area_km2 = target_area,
    n_uci_holdout_in_area_band = nrow(matched),
    n_cta_queries = nrow(target),
    uci_baseline_mae_log1p = if (nrow(matched) >= 10) {
      mean(matched$abs_log_error_baseline)
    } else NA_real_,
    uci_augmented_mae_log1p = if (nrow(matched) >= 10) {
      mean(matched$abs_log_error_augmented)
    } else NA_real_,
    cta_baseline_mae_log1p = mean(target$abs_log_error_baseline),
    cta_augmented_mae_log1p = mean(target$abs_log_error_augmented),
    stringsAsFactors = FALSE
  )
}))
write.csv(area_checks, file.path(output_dir, '05_area_overlap_diagnostic.csv'),
          row.names = FALSE)

# Geographic/area support diagnostic: how many CTA queries are in the
# training set's central 98% range for area AND both center coordinates?
within_range <- function(v, reference) {
  q <- quantile(reference, c(0.01, 0.99), na.rm = TRUE)
  v >= q[1] & v <= q[2]
}
cta$within_uci_training_support <-
  within_range(cta$log_area, train$log_area) &
  within_range(cta$center_lat, train$center_lat) &
  within_range(cta$center_lon, train$center_lon)

# Select a preferred model based ONLY on the UCI held-out MAE in log space.
holdout_scores <- evaluation[evaluation$evaluation_set == 'UCI_spatial_holdout', ]
preferred <- holdout_scores$model[which.min(holdout_scores$mae_log1p_count)]

# Paired difference baseline - augmented (positive means augmentation helped).
# Bootstrap samples CTA STATIONS, not individual queries, because each station
# contributes three correlated radii. Nearby stations can still overlap.
set.seed(5712025)
cta$mae_log_gain <- cta$abs_log_error_baseline - cta$abs_log_error_augmented
station_values <- split(cta$mae_log_gain, as.character(cta$station_id))
station_ids <- names(station_values)
if (length(station_ids) < 20L) stop('Too few CTA stations for bootstrap.')
bootstrap_mean <- replicate(1000L, {
  ids <- sample(station_ids, length(station_ids), replace = TRUE)
  mean(unlist(station_values[ids], use.names = FALSE))
})
boot_ci <- as.numeric(quantile(bootstrap_mean, c(0.025, 0.975)))
mean_gain <- mean(cta$mae_log_gain)

# --- 4. Save reproducible predictions and model parameters -------------
keep <- c('query_id', 'query_source', 'station_id', 'query_halfwidth_m',
          'area_km2', 'center_lat', 'center_lon', 'stations_inside_query',
          'total_cta_entries_at_stations_inside_query_2024',
          'crimes_2024', 'pred_baseline_count', 'pred_augmented_count',
          'abs_log_error_baseline', 'abs_log_error_augmented')
predictions <- rbind(holdout[, keep], cta[, keep])
predictions$evaluation_set <- c(rep('UCI_spatial_holdout', nrow(holdout)),
                                rep('CTA_transfer', nrow(cta)))
predictions$within_uci_training_support <-
  c(rep(NA, nrow(holdout)), cta$within_uci_training_support)
write.csv(predictions,
          file.path(processed_dir, '05_model_predictions_2024.csv'),
          row.names = FALSE)
coefs <- summary(augmented)$coefficients
write.csv(data.frame(term = rownames(coefs), coefs, row.names = NULL,
                     check.names = FALSE),
          file.path(output_dir, '05_augmented_model_coefficients.csv'),
          row.names = FALSE)

# --- 5. Three plots ----------------------------------------------------
# Chart 1: compare log-scale MAE across UCI holdout and each CTA size.
plot_sets <- c('UCI_spatial_holdout', 'CTA_250m_square',
               'CTA_500m_square', 'CTA_750m_square')
metric_mat <- sapply(plot_sets, function(s) {
  sub <- evaluation[evaluation$evaluation_set == s, ]
  sub$mae_log1p_count[match(c('Baseline', 'Augmented'), sub$model)]
})
rownames(metric_mat) <- c('Area + geography', 'Area + geography + CTA')
colnames(metric_mat) <- c('UCI holdout', 'CTA 250m', 'CTA 500m', 'CTA 750m')
png(file.path(output_dir, '05_transfer_error_comparison.png'),
    width = 1360, height = 760, res = 125)
op <- par(mar = c(6, 5, 4, 2))
barplot(metric_mat, beside = TRUE, col = c('grey40', 'grey75'),
        ylim = c(0, max(metric_mat) * 1.18), las = 1,
        main = 'Generalization: UCI spatial holdout vs CTA-centered queries',
        ylab = 'Mean absolute error of log(1 + 2024 crime count)',
        legend.text = rownames(metric_mat),
        args.legend = list(x = 'topright', bty = 'n', cex = 0.8))
par(op)
dev.off()

# Chart 2: observed versus estimated query answers on CTA transfer set.
png(file.path(output_dir, '05_cta_predicted_vs_observed.png'),
    width = 1400, height = 740, res = 120)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
limit <- range(c(log1p(cta$crimes_2024), cta$pred_baseline_log,
                 cta$pred_augmented_log), finite = TRUE)
for (method in c('Baseline', 'Augmented')) {
  pred <- if (method == 'Baseline') cta$pred_baseline_log else cta$pred_augmented_log
  plot(log1p(cta$crimes_2024), pred, pch = 19, cex = 0.65,
       col = 'grey45', xlim = limit, ylim = limit,
       xlab = 'Observed log(1 + 2024 crime count)',
       ylab = 'Predicted log(1 + 2024 crime count)',
       main = paste(method, 'model: CTA transfer'))
  abline(0, 1, lty = 2, col = 'grey20', lwd = 1.5)
}
par(op)
dev.off()

# Chart 3: where CTA feature augmentation helps or hurts by query size.
png(file.path(output_dir, '05_cta_feature_gain_by_radius.png'),
    width = 1000, height = 720, res = 120)
op <- par(mar = c(5, 5, 4, 2))
gain <- sapply(c(250, 500, 750), function(radius) {
  mean(cta$mae_log_gain[cta$query_halfwidth_m == radius])
})
barplot(gain, names.arg = paste(c(250, 500, 750), 'm'),
        col = 'grey60', border = NA,
        main = 'CTA features: change in prediction error by square halfwidth',
        xlab = 'CTA-centered square halfwidth',
        ylab = 'Baseline MAE minus augmented MAE (log scale)')
abline(h = 0, lty = 2)
par(op)
dev.off()

# --- 6. Summary: don't claim causality or passenger risk ---------------
report <- c(
  'CSP 571 - THE TRANSIT EXPOSURE PARADOX',
  'STEP 5: QUERY ESTIMATION AND SYNTHETIC-TO-TRANSIT TRANSFER',
  paste('Generated:', Sys.time()),
  paste('Training UCI queries:', nrow(train)),
  paste('Geographically held-out UCI queries:', nrow(holdout)),
  paste('CTA transfer queries:', nrow(cta),
        '(across', length(station_ids), 'station IDs and three square sizes)'),
  paste('CTA queries within central 98% UCI training support (area + lon + lat):',
        sum(cta$within_uci_training_support), 'of', nrow(cta)),
  '',
  'MODELS (both fitted ONLY on UCI training rectangles):',
  paste('Baseline:', deparse(baseline_formula)),
  paste('Augmented:', paste(deparse(augmented_formula), collapse = ' ')),
  paste('Preferred by UCI geographically held-out MAE_log1p:', preferred),
  '',
  'PREDICTION METRICS:',
  capture.output(print(evaluation, row.names = FALSE, digits = 4)),
  '',
  'AREA-OVERLAP DIAGNOSTIC (area within factor 1.5 of CTA square):',
  capture.output(print(area_checks, row.names = FALSE, digits = 4)),
  '',
  'CTA PAIRED FEATURE EFFECT:',
  paste('Mean baseline-minus-augmented absolute log error:', signif(mean_gain, 4)),
  paste('Approximate station-cluster bootstrap 95% CI:',
        paste(signif(boot_ci, 4), collapse = ' to ')),
  '(Positive means the augmented model makes smaller errors.)',
  '',
  'CRITICAL INTERPRETATION:',
  '- Outcome is exact 2024 spatial-query count, not historical UCI count.',
  '- The 800 UCI queries were sampled partly by area; they are not a random sample of all UCI records.',
  '- UCI test rectangles hold out center-location tiles, but polygons can still overlap tile borders.',
  '- CTA samples were NEVER used to fit either model; this is an out-of-workload test.',
  '- Location distributions, query areas, and rectangle aspect ratios differ across workloads.',
  '- Log-scale errors permit comparisons across count scales; matching area alone does not control geography.',
  '- Current-year CTA ridership measures activity, not personal crime exposure or causal crime risk.',
  '- This is contemporaneous query estimation, NOT forecasting future crimes.',
  '- Bootstrap clusters by station, but neighboring station queries may still overlap geographically.',
  '- A smaller held-out error suggests better estimation only for the tested sample and split.',
  '',
  'OUTPUT FILES:',
  'data/processed/05_model_predictions_2024.csv',
  'outputs/05_uci_spatial_split.csv',
  'outputs/05_evaluation_metrics.csv',
  'outputs/05_area_overlap_diagnostic.csv',
  'outputs/05_augmented_model_coefficients.csv',
  'outputs/05_transfer_error_comparison.png',
  'outputs/05_cta_predicted_vs_observed.png',
  'outputs/05_cta_feature_gain_by_radius.png'
)
writeLines(report, file.path(output_dir, '05_benchmark_transfer_summary_2024.txt'))

cat('\nSUCCESS - STEP 5 benchmark transfer completed.\n')
cat('Preferred UCI model: ', preferred, '\n', sep = '')
cat('CTA paired effect (positive is better): ', signif(mean_gain, 4), '\n', sep = '')
cat('CTA bootstrap interval: ', paste(signif(boot_ci, 4), collapse = ' to '), '\n', sep = '')
cat('Results saved in outputs/ and data/processed/.\n')
