# CSP 571 - Transit Exposure Paradox
# Step 7A: audit UCI coverage, integrate Task A, preserve Step 5 split.
# Run from repository root. Requires sf. No model fitting in this step.
if (!requireNamespace('sf', quietly = TRUE)) stop('Install sf first: install.packages("sf")')
library(sf)
sf_use_s2(FALSE)
options(stringsAsFactors = FALSE)
read_input <- function(path) {
  if (!file.exists(path)) stop('Missing input: ', path, '\nRun from the repository root.')
  read.csv(path, stringsAsFactors = FALSE)
}
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
unique_ids <- function(x, label) assert(!anyNA(x) && !anyDuplicated(x), paste('Invalid or duplicate IDs:', label))
read_num <- function(x) suppressWarnings(as.numeric(x))

w <- read_input('data/processed/04b_spatial_query_workloads_2024.csv')
s <- read_input('outputs/04b_uci_query_sample.csv')
sp <- read_input('outputs/05_uci_spatial_split.csv')
cov <- read_input('data/processed/cta_geographic_coverage.csv')
pred <- read_input('data/processed/05_model_predictions_2024.csv')
unique_ids(w$query_id, 'workload'); unique_ids(sp$query_id, 'split')
unique_ids(cov$query_id, 'CTA coverage'); unique_ids(pred$query_id, 'predictions')
uci <- w[w$query_source == 'UCI_original_rectangle', , drop = FALSE]
cta <- w[w$query_source == 'CTA_centered_rectangle', , drop = FALSE]
assert(nrow(uci) == 800L && nrow(cta) == 426L && nrow(w) == 1226L, 'Unexpected Step 4B workload sizes.')
s$query_id <- paste0('UCI_ROW_', s$source_row_number)
unique_ids(s$query_id, 'UCI sample')
assert(setequal(s$query_id, uci$query_id) && setequal(sp$query_id, uci$query_id), 'UCI sample/split IDs differ from workload.')
assert(setequal(cov$query_id, cta$query_id), 'Task A coverage IDs differ from CTA workload.')
s <- s[match(uci$query_id, s$query_id), , drop = FALSE]
sp <- sp[match(uci$query_id, sp$query_id), , drop = FALSE]
cov <- cov[match(cta$query_id, cov$query_id), , drop = FALSE]
assert(all(sp$split %in% c('UCI_train', 'UCI_spatial_holdout')), 'Unknown Step 5 split label.')
assert(sum(sp$split == 'UCI_train') == 660L && sum(sp$split == 'UCI_spatial_holdout') == 140L, 'Unexpected Step 5 split counts.')
assert(!anyNA(sp$spatial_block) && !length(intersect(sp$spatial_block[sp$split == 'UCI_train'], sp$spatial_block[sp$split == 'UCI_spatial_holdout'])), 'Spatial block crosses split.')
assert(all(cov$crimes_2024 == cta$crimes_2024), 'CTA coverage counts do not match workload.')

# Same UCI reconstruction as Step 4B. EPSG:3435 is the inferred CRS used there;
# this step does not independently establish the original benchmark CRS.
rect <- function(x, y, dx, dy) {
  st_polygon(list(matrix(c(x-dx,y-dy, x+dx,y-dy, x+dx,y+dy,
                           x-dx,y+dy, x-dx,y-dy), ncol = 2, byrow = TRUE)))
}
cols <- c('center_x_stateplane_ft','center_y_stateplane_ft','halfwidth_x_ft','halfheight_y_ft')
assert(all(cols %in% names(s)), 'Missing original rectangle coordinate fields.')
for (nm in cols) s[[nm]] <- read_num(s[[nm]])
assert(all(is.finite(as.matrix(s[,cols]))) && all(s$halfwidth_x_ft > 0 & s$halfheight_y_ft > 0), 'Invalid UCI coordinates or ranges.')
uci_geom <- st_transform(st_sfc(lapply(seq_len(nrow(s)), function(i) {
  rect(s$center_x_stateplane_ft[i], s$center_y_stateplane_ft[i], s$halfwidth_x_ft[i], s$halfheight_y_ft[i])
}), crs = 3435), 26916)
uci_area <- as.numeric(st_area(uci_geom))
assert(all(abs(uci_area / 1e6 - uci$area_km2) <= pmax(1e-8, 1e-6 * uci$area_km2)), 'Reconstructed UCI areas differ from Step 4B.')

# Reuse the exact committed Task A boundary snapshot; avoid local raw-file ambiguity.
boundary_path <- 'data/Boundaries_-_City_20261121.geojson'
assert(file.exists(boundary_path), paste('Missing committed boundary:', boundary_path))
city <- st_read(boundary_path, quiet = TRUE)
city_metric <- st_make_valid(st_transform(st_union(st_make_valid(city)), 26916))
assert(length(city_metric) == 1L && all(st_is_valid(city_metric)) && !any(st_is_empty(city_metric)), 'Invalid Chicago boundary.')
cat('Intersecting 800 UCI rectangle footprints with Chicago city boundary...\n')
q_sf <- st_sf(row_idx = seq_len(nrow(uci)), geometry = uci_geom)
overlap <- suppressWarnings(st_intersection(q_sf, city_metric))
inside <- numeric(nrow(uci))
if (nrow(overlap) > 0) {
  sums <- tapply(as.numeric(st_area(overlap)), overlap$row_idx, sum)
  inside[as.integer(names(sums))] <- as.numeric(sums)
}
assert(all(is.finite(inside)) && all(inside >= 0) && all(inside <= uci_area + pmax(1e-4, 1e-8 * uci_area)), 'Intersection areas invalid.')
# Use Task A's operational tolerances (percent, not fraction).
classify <- function(p) ifelse(p >= 99.99, 'fully covered', ifelse(p <= 0.0001, 'outside Chicago', 'partially covered'))
uci_pct <- pmin(100, pmax(0, 100 * inside / uci_area))
cta_pct <- read_num(cov$pct_inside_chicago)
assert(all(is.finite(cta_pct) & cta_pct >= 0 & cta_pct <= 100), 'Invalid CTA coverage percentage.')
assert(all(classify(cta_pct) == cov$coverage_status), 'CTA coverage labels inconsistent with Task A thresholds.')
assert(sum(cov$coverage_status == 'fully covered') == 340L && sum(cov$coverage_status == 'partially covered') == 36L && sum(cov$coverage_status == 'outside Chicago') == 50L, 'Unexpected Task A classification totals.')
assert(sum(cta$crimes_2024 == 0) == 54L && !any(cta$crimes_2024 == 0 & cov$coverage_status == 'fully covered'), 'Unexpected CTA zero-count results.')

# One row per original query; retain all rows and all zeros, even outside city.
# Coverage is a jurisdictional proxy, NOT a guarantee of complete incident data.
uci$split <- sp$split; uci$spatial_block <- sp$spatial_block
uci$pct_inside_chicago <- uci_pct; uci$coverage_status <- classify(uci_pct)
cta$split <- 'CTA_transfer'; cta$spatial_block <- NA_character_
cta$pct_inside_chicago <- cta_pct; cta$coverage_status <- cov$coverage_status
all <- rbind(uci, cta)
all$inside_city_area_km2 <- c(inside / 1e6, cov$inside_area_m2 / 1e6)
all$zero_crime_query <- all$crimes_2024 == 0
all$eligible_full_coverage <- all$coverage_status == 'fully covered'
unique_ids(all$query_id, 'integrated workload')
assert(setequal(all$query_id, w$query_id) && nrow(all) == 1226L, 'Integrated workload lost queries.')
assert(!anyNA(all$coverage_status) && !anyNA(all$zero_crime_query), 'Missing coverage or zero flag.')
all <- all[match(w$query_id, all$query_id), , drop = FALSE]

# Coverage summary: no count-based filtering. Zeros in fully covered areas stay eligible.
classes <- c('fully covered','partially covered','outside Chicago')
splits <- c('UCI_train','UCI_spatial_holdout','CTA_transfer')
summary <- do.call(rbind, lapply(splits, function(g) {
  a <- all[all$split == g, , drop = FALSE]
  do.call(rbind, lapply(classes, function(k) {
    z <- a[a$coverage_status == k, , drop = FALSE]
    data.frame(split = g, coverage_status = k, n_queries = nrow(z),
               pct_of_split = 100 * nrow(z)/nrow(a), zero_queries = sum(z$zero_crime_query),
               median_area_km2 = if (nrow(z)) median(z$area_km2) else NA_real_)
  }))
}))

# Evaluate EXISTING Step 5 predictions, no retraining yet. This isolates
# differences across coverage strata for a fixed fitted model; not causal.
expected <- all$query_id[all$split != 'UCI_train']
assert(setequal(pred$query_id, expected), 'Prediction IDs differ from Step 5 test split.')
m <- match(pred$query_id, all$query_id)
assert(all(pred$crimes_2024 == all$crimes_2024[m]) && all(pred$evaluation_set == all$split[m]), 'Prediction alignment failed.')
pred$coverage_status <- all$coverage_status[m]
pred$pct_inside_chicago <- all$pct_inside_chicago[m]
metrics <- do.call(rbind, lapply(splits[-1], function(g) {
  a <- pred[pred$evaluation_set == g, , drop = FALSE]
  do.call(rbind, lapply(c('all queries', classes), function(k) {
    z <- if (k == 'all queries') a else a[a$coverage_status == k, , drop = FALSE]
    do.call(rbind, lapply(c('baseline','augmented'), function(model) {
      e <- z[[paste0('abs_log_error_', model)]]
      p <- z[[paste0('pred_', model, '_count')]]
      assert(all(is.finite(e)) && all(is.finite(p)), 'Nonfinite prediction values.')
      data.frame(evaluation_set = g, coverage_status = k, model = model,
                 n_queries = nrow(z), mae_log1p = if (nrow(z)) mean(e) else NA_real_,
                 mae_count = if (nrow(z)) mean(abs(p - z$crimes_2024)) else NA_real_)
    }))
  }))
}))

out <- 'outputs'; dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(all, 'data/processed/07a_coverage_aware_workloads_2024.csv', row.names = FALSE)
write.csv(summary, file.path(out, '07a_coverage_summary.csv'), row.names = FALSE)
write.csv(metrics, file.path(out, '07a_existing_model_coverage_metrics.csv'), row.names = FALSE)
full_train <- all[all$split == 'UCI_train' & all$eligible_full_coverage, , drop = FALSE]
full_holdout <- all[all$split == 'UCI_spatial_holdout' & all$eligible_full_coverage, , drop = FALSE]
full_cta <- all[all$split == 'CTA_transfer' & all$eligible_full_coverage, , drop = FALSE]
block_n <- length(unique(full_train$spatial_block))
cat('\nCOVERAGE BY ORIGINAL SPLIT\n'); print(summary, row.names = FALSE, digits = 4)
cat('\nEXISTING MODEL ERRORS BY COVERAGE\n'); print(metrics, row.names = FALSE, digits = 4)
lines <- c('STEP 7A - COVERAGE-AWARE WORKLOAD INTEGRATION',
  paste('Boundary:', boundary_path),
  'UCI CRS assumed EPSG:3435 as in Step 4B; intersections use EPSG:26916.',
  'Coverage tolerance: >=99.99% fully covered; <=0.0001% outside; otherwise partial.',
  'All original queries, counts and Step 5 spatial split retained.',
  '', 'COVERAGE BY SPLIT:', capture.output(print(summary, row.names = FALSE, digits = 4)),
  '', 'EXISTING STEP 5 MODEL ERRORS (no refitting):', capture.output(print(metrics, row.names = FALSE, digits = 4)),
  '', 'FULL-COVERAGE MODELING CANDIDATES:',
  paste('UCI training:', nrow(full_train), 'queries across', block_n, 'center blocks'),
  paste('UCI holdout:', nrow(full_holdout)), paste('CTA transfer:', nrow(full_cta)),
  'Use only UCI training rows for tuning/fitting. Keep the original holdout and CTA transfer untouched.',
  'Check retained training size and block count before choosing CV folds or nonlinear model complexity.',
  'Municipal boundary is a coverage proxy; it does not prove complete reporting.',
  'UCI center-block splitting does not guarantee nonoverlapping query polygons.',
  'Differences across coverage groups are descriptive, not causal. No zeros were discarded based on count.',
  'No model was fitted or selected in Step 7A.')
writeLines(lines, file.path(out, '07a_integration_summary.txt'))
png(file.path(out, '07a_coverage_by_workload.png'), width = 1250, height = 720, res = 120)
mat <- xtabs(n_queries ~ coverage_status + split, data = summary)
mat <- mat[classes, splits, drop = FALSE]
par(mar = c(5,5,4,2))
barplot(prop.table(mat, 2)*100, col = c('#238b8d','#f2b447','#ba5656'), border = NA,
        names.arg = c('UCI train','UCI holdout','CTA transfer'), ylim = c(0, 120),
        ylab = 'Percent of original split', main = 'Chicago geographic coverage by workload',
        legend.text = c('Fully covered','Partially covered','Outside Chicago'),
        args.legend = list(x = 'top', horiz = TRUE, bty = 'n', cex = 0.85))
dev.off()
cat('\nSTEP 7A COMPLETE\n')
cat('Fully covered UCI train:', nrow(full_train), 'across', block_n, 'blocks\n')
cat('Fully covered UCI holdout:', nrow(full_holdout), '\n')
cat('Fully covered CTA transfer:', nrow(full_cta), '\n')
cat('Review outputs/07a_integration_summary.txt before model training.\n')
