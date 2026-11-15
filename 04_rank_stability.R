# CSP 571 - The Transit Exposure Paradox
# STEP 4A: Rank reversals and sensitivity to spatial query radius
# Input: data/processed/station_crime_ridership_2024.csv (from Step 3)
# Base R only; no new downloads, external packages, or spatial joins.

input_path <- file.path('data', 'processed', 'station_crime_ridership_2024.csv')
output_dir <- 'outputs'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop('Missing input: ', input_path,
    '\nRun this script from the project root. Current directory: ', getwd())

x <- read.csv(input_path, stringsAsFactors = FALSE)
needed <- c('station_id','station_name','query_radius_m',
            'total_entries_2024','reported_crimes_2024',
            'crimes_per_100k_entries','raw_count_rank','activity_adjusted_rank')
missing <- setdiff(needed, names(x))
if (length(missing)) stop('Missing columns: ', paste(missing, collapse = ', '))

# Each station should have exactly one row per radius. Group by ID, never name
# (separate stations can share names, such as Jackson and Monroe).
x$station_id <- as.character(x$station_id)
x$query_radius_m <- as.integer(x$query_radius_m)
x$station_name <- as.character(x$station_name)
for (nm in c('total_entries_2024','reported_crimes_2024',
             'crimes_per_100k_entries','raw_count_rank','activity_adjusted_rank')) {
  x[[nm]] <- as.numeric(x[[nm]])
}
radii <- c(250L, 500L, 750L)
if (!identical(sort(unique(x$query_radius_m)), radii)) {
  stop('Expected exactly these query radii: 250, 500, 750 meters.')
}
key <- paste(x$station_id, x$query_radius_m, sep = '|')
if (anyDuplicated(key)) stop('Duplicate station ID + radius detected.')
if (anyNA(x[, needed]) || any(x$total_entries_2024 <= 0) ||
    any(x$reported_crimes_2024 < 0)) {
  stop('Missing or invalid required fields; inspect Step 3 output.')
}
ids_by_radius <- lapply(radii, function(r) sort(x$station_id[x$query_radius_m == r]))
if (!all(vapply(ids_by_radius[-1], identical, logical(1), ids_by_radius[[1]]))) {
  stop('Station ID sets differ between query radii; do not compare unmatched ranks.')
}

# A station moving from raw rank 1 to adjusted rank 40 has a rank shift of +39:
# a positive shift means it moves DOWN the highest-count ranking after adjustment.
x$adjusted_minus_raw_rank <- x$activity_adjusted_rank - x$raw_count_rank
x$absolute_rank_shift <- abs(x$adjusted_minus_raw_rank)
x <- x[order(x$query_radius_m, x$raw_count_rank, x$station_id), ]
write.csv(x, file.path(output_dir, '04_station_rank_changes_all_radii.csv'),
          row.names = FALSE)

# Within each radius, compare raw-crime-count rankings to ridership-adjusted rankings.
within <- do.call(rbind, lapply(radii, function(r) {
  d <- x[x$query_radius_m == r, ]
  data.frame(
    radius_m = r,
    stations = nrow(d),
    spearman_raw_vs_adjusted = cor(d$raw_count_rank,
                                   d$activity_adjusted_rank, method = 'spearman'),
    median_absolute_rank_shift = median(d$absolute_rank_shift),
    mean_absolute_rank_shift = mean(d$absolute_rank_shift),
    stations_shifted_20_or_more_places = sum(d$absolute_rank_shift >= 20),
    zero_crime_station_queries = sum(d$reported_crimes_2024 == 0)
  )
}))
write.csv(within, file.path(output_dir, '04_within_radius_metrics.csv'),
          row.names = FALSE)

# Across radii, compare the same stations by their ranks. Counts increase as the
# buffer grows, so ranking comparisons are more interpretable than raw counts.
get_r <- function(r) {
  d <- x[x$query_radius_m == r, ]
  d[order(d$station_id), ]
}
top_n <- 10L
get_top <- function(d, metric) {
  d$station_id[head(order(d[[metric]], d$station_id), top_n)]
}
rad_pairs <- combn(radii, 2L, simplify = FALSE)
across <- do.call(rbind, lapply(rad_pairs, function(p) {
  a <- get_r(p[1]); b <- get_r(p[2])
  stopifnot(identical(a$station_id, b$station_id))
  top_raw_a <- get_top(a, 'raw_count_rank')
  top_raw_b <- get_top(b, 'raw_count_rank')
  top_adj_a <- get_top(a, 'activity_adjusted_rank')
  top_adj_b <- get_top(b, 'activity_adjusted_rank')
  overlap_raw <- length(intersect(top_raw_a, top_raw_b))
  overlap_adj <- length(intersect(top_adj_a, top_adj_b))
  data.frame(
    radius_a_m = p[1], radius_b_m = p[2],
    spearman_raw_count_ranks = cor(a$raw_count_rank, b$raw_count_rank,
                                   method = 'spearman'),
    spearman_adjusted_ranks = cor(a$activity_adjusted_rank,
                                  b$activity_adjusted_rank, method = 'spearman'),
    raw_top10_overlap = overlap_raw,
    adjusted_top10_overlap = overlap_adj,
    raw_top10_jaccard = overlap_raw / (2 * top_n - overlap_raw),
    adjusted_top10_jaccard = overlap_adj / (2 * top_n - overlap_adj)
  )
}))
write.csv(across, file.path(output_dir, '04_cross_radius_metrics.csv'),
          row.names = FALSE)

# A quick table of the most striking 500 m changes, with IDs to resolve
# repeated station names.
d500 <- get_r(500L)
d500$station_label <- paste0(d500$station_name, ' [', d500$station_id, ']')
movers <- d500[order(-d500$absolute_rank_shift, d500$station_id), ]
write.csv(movers[, c('station_id','station_name',
                     'reported_crimes_2024','total_entries_2024',
                     'crimes_per_100k_entries','raw_count_rank',
                     'activity_adjusted_rank','adjusted_minus_raw_rank')],
          file.path(output_dir, '04_largest_rank_shifts_500m.csv'),
          row.names = FALSE)

# ---- Figure 1: 500m raw versus adjusted rankings ----------------------
png(file.path(output_dir, '04_raw_vs_adjusted_ranks_500m.png'),
    width = 1150, height = 950, res = 125)
op <- par(mar = c(5, 5, 4, 2))
plot(d500$raw_count_rank, d500$activity_adjusted_rank,
     pch = 19, cex = 0.65, col = 'grey40',
     xlab = 'Raw crime-count rank (1 = highest count)',
     ylab = 'Ridership-adjusted rank (1 = highest ratio)',
     main = 'CTA station rank changes after ridership adjustment: 500 m',
     xlim = c(1, nrow(d500)), ylim = c(nrow(d500), 1))
abline(a = 0, b = 1, lty = 2, col = 'grey60', lwd = 2)
legend('bottomleft', legend = 'Dashed line = unchanged rank', lty = 2,
       col = 'grey60', bty = 'n', cex = 0.85)
par(op)
dev.off()

# ---- Figure 2: cross-radius ranking stability ------------------------
a <- get_r(250L); b <- get_r(750L)
png(file.path(output_dir, '04_250m_vs_750m_rank_stability.png'),
    width = 1500, height = 750, res = 120)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
for (metric in c('raw_count_rank', 'activity_adjusted_rank')) {
  plot(a[[metric]], b[[metric]], pch = 19, cex = 0.62,
       col = 'grey45',
       xlab = '250 m rank (1 = highest)',
       ylab = '750 m rank (1 = highest)',
       main = if (metric == 'raw_count_rank') 'Raw crime-count rankings'
              else 'Ridership-adjusted rankings',
       xlim = c(1, nrow(a)), ylim = c(nrow(a), 1))
  abline(a = 0, b = 1, lty = 2, col = 'grey60', lwd = 2)
}
par(op)
dev.off()

# ---- Figure 3: biggest 500m rank movers -------------------------------
show <- head(movers, 12L)
png(file.path(output_dir, '04_largest_rank_shifts_500m.png'),
    width = 1320, height = 850, res = 120)
op <- par(mar = c(5, 19, 4, 2))
barplot(rev(show$adjusted_minus_raw_rank),
        names.arg = rev(show$station_label), horiz = TRUE, las = 1,
        col = 'grey55', border = NA,
        xlab = 'Adjusted rank minus raw rank (+ = moves down the list)',
        main = 'Largest changes in station rankings: 500 m')
abline(v = 0, lty = 2, col = 'grey35')
par(op)
dev.off()

# ---- Text summary, including non-causal interpretation ----------------
report <- c(
  'CSP 571 - THE TRANSIT EXPOSURE PARADOX',
  'STEP 4A: SPATIAL QUERY RANKING SENSITIVITY',
  paste('Generated:', Sys.time()),
  paste('Stations compared at all radii:', nrow(d500)),
  'Radii: 250 m, 500 m, 750 m', '',
  'WITHIN-RADIUS: RAW VS RIDERSHIP-ADJUSTED RANKS',
  capture.output(print(within, row.names = FALSE, digits = 4)), '',
  'CROSS-RADIUS: RANK CONSISTENCY AND TOP-10 OVERLAP',
  capture.output(print(across, row.names = FALSE, digits = 4)), '',
  'TOP 12 ABSOLUTE RANK SHIFTS AT 500 M',
  capture.output(print(head(movers[, c('station_label', 'raw_count_rank',
                             'activity_adjusted_rank', 'adjusted_minus_raw_rank')], 12),
                       row.names = FALSE)), '',
  'IMPORTANT INTERPRETATION NOTES',
  '- This measures changes in query rankings; it does NOT estimate an individual passenger\'s crime risk.',
  '- CTA entries are an imperfect indicator of nearby pedestrian activity.',
  '- Nearby reported incidents are not necessarily transit-related.',
  '- A change of ranking after dividing by ridership is not proof that either ranking is correct.',
  '- Some station buffers overlap and some stations lie outside the Chicago crime-data coverage.',
  '- Crime geographic coordinates are approximate; 250m query counts are particularly sensitive.',
  '- Results are descriptive for 2024, not causal or predictive.',
  '- Original UCI spatial queries have not yet been linked; that is Step 4B.'
)
writeLines(report, file.path(output_dir, '04_rank_stability_summary_2024.txt'))

cat('\nSUCCESS - Step 4A ranking sensitivity complete.\n')
cat('Stations compared:', nrow(d500), '\n')
cat('Saved 4 CSVs, 3 plots and 1 summary text in outputs/.\n')
print(within, row.names = FALSE)
print(across, row.names = FALSE)
