# CSP 571 - The Transit Exposure Paradox
# STEP 4B: Augment original UCI rectangular spatial queries with 2024
# Chicago crime and CTA station ridership; construct comparable CTA-centered
# RECTANGULAR queries for workload-distribution analysis.
#
# Run from your project root AFTER scripts 01, 02, and 03.
# Requirements: install.packages("sf") (if not installed already).
# No APIs or downloads required.
#
# IMPORTANT:
# - UCI's Range-Queries-Aggregates.csv uses x +/- x_range and y +/- y_range.
# - The X/Y values are consistent with Chicago's Illinois East State Plane
#   coordinates (NAD83, US survey feet; EPSG:3435). We sanity-check this.
# - UCI's two radius-query files use MUCH SMALLER undocumented coordinate
#   values. DO NOT overlay them on Chicago without an independently validated
#   transformation. This script intentionally uses UCI's RECTANGULAR file.
# - Old UCI aggregate counts and 2024 counts are different historical
#   snapshots. They are preserved as separate fields, NOT compared as errors.

if (!requireNamespace("sf", quietly = TRUE)) {
  stop('Package "sf" is required. Run install.packages("sf") once.')
}
library(sf)

out_dir <- "outputs"
processed_dir <- file.path("data", "processed")
raw_uci_dir <- file.path("data", "raw", "uci")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

crime_path <- file.path(processed_dir, "chicago_crimes_2024_clean.csv")
station_path <- file.path(processed_dir, "station_crime_ridership_2024.csv")
uci_candidates <- list.files(raw_uci_dir, pattern = "^Range-Queries-Aggregates\\.csv$",
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = TRUE)
if (length(uci_candidates) != 1L) {
  stop("Expected one Range-Queries-Aggregates.csv under data/raw/uci/. ",
       "Rerun 01_inspect_uci.R if it is missing.")
}
for (p in c(crime_path, station_path)) {
  if (!file.exists(p)) stop("Missing input file: ", p)
}

cat("Reading original UCI rectangular queries and processed 2024 files...\n")
uci <- read.csv(uci_candidates[1], check.names = FALSE,
                stringsAsFactors = FALSE)
crime <- read.csv(crime_path, stringsAsFactors = FALSE)
all_stations <- read.csv(station_path, stringsAsFactors = FALSE)

# Accept the published UCI file's x, y, x_range, y_range, count, sum_, avg
# (some versions include an unnamed leading row-ID column).
names(uci) <- tolower(trimws(names(uci)))
required_uci <- c("x", "y", "x_range", "y_range", "count", "sum_", "avg")
if (!all(required_uci %in% names(uci))) {
  stop("Unexpected UCI columns: ", paste(names(uci), collapse = ", "))
}
required_crime <- c("latitude", "longitude", "primary_type")
required_station <- c("station_id", "station_name", "latitude", "longitude",
                      "total_entries_2024", "query_radius_m")
if (!all(required_crime %in% names(crime))) stop("Processed crime schema changed.")
if (!all(required_station %in% names(all_stations))) {
  stop("Processed station schema changed. Rerun Step 3.")
}

# Use the same 142 station IDs from Step 3, one row per physical station.
stations <- all_stations[all_stations$query_radius_m == 500, ]
if (!nrow(stations) || anyDuplicated(stations$station_id)) {
  stop("Station IDs missing/duplicated at 500m; inspect Step 3 output.")
}
if (!nrow(crime)) stop("No cleaned crime records.")

# Reproducible UCI sample: 600 arbitrary queries + up to 200 roughly
# area-matched queries (near the area of a 500m-halfwidth CTA square).
# Sampling is used for feasible local spatial processing, NOT for a
# representative claim about every 200,000 UCI queries.
for (v in c("x", "y", "x_range", "y_range", "count", "sum_", "avg")) {
  uci[[v]] <- suppressWarnings(as.numeric(uci[[v]]))
}
valid <- is.finite(uci$x) & is.finite(uci$y) &
  is.finite(uci$x_range) & is.finite(uci$y_range) &
  uci$x_range > 0 & uci$y_range > 0
uci <- uci[valid, ]
uci$source_row_number <- which(valid)
if (nrow(uci) < 100) stop("Too few valid UCI queries.")

# In State Plane Illinois East EPSG:3435, distance is in US survey feet.
ft_to_m <- 1200 / 3937
uci$approx_area_km2 <- 4 * uci$x_range * uci$y_range * ft_to_m^2 / 1e6
set.seed(5712024)
random_rows <- sample.int(nrow(uci), min(600L, nrow(uci)))
size_candidates <- which(uci$approx_area_km2 >= 0.25 &
                         uci$approx_area_km2 <= 4.0)
size_candidates <- setdiff(size_candidates, random_rows)
size_rows <- if (length(size_candidates)) {
  sample(size_candidates, min(200L, length(size_candidates)))
} else integer(0)
selected <- unique(c(random_rows, size_rows))
u <- uci[selected, , drop = FALSE]
rownames(u) <- NULL
cat("Sampled UCI queries:", nrow(u), "(random:", length(random_rows),
    "; area-targeted:", length(size_rows), ")\n")

# Verify plausible geographic location BEFORE spatial intersections.
centers <- st_as_sf(u, coords = c("x", "y"), crs = 3435,
                    remove = FALSE)
lonlat <- st_coordinates(st_transform(centers, 4326))
plausible <- is.finite(lonlat[, 1]) & is.finite(lonlat[, 2]) &
             lonlat[, 1] > -88.2 & lonlat[, 1] < -87.2 &
             lonlat[, 2] > 41.5 & lonlat[, 2] < 42.2
cat("UCI query centers inside broad Chicago-area check:",
    sum(plausible), "/", length(plausible), "\n")
if (mean(plausible) < 0.80) {
  stop("Coordinate validation FAILED; do not assume UCI EPSG:3435.")
}

# Rectangle maker: UCI explicitly defines x +/- x_range, y +/- y_range.
make_rectangle <- function(xmin, ymin, xmax, ymax) {
  st_polygon(list(matrix(c(xmin, ymin,
                           xmax, ymin,
                           xmax, ymax,
                           xmin, ymax,
                           xmin, ymin),
                         ncol = 2, byrow = TRUE)))
}

# Original UCI rectangles are in survey feet. Convert polygon to UTM meters.
uci_rect_ft <- st_sfc(lapply(seq_len(nrow(u)), function(i) {
  make_rectangle(u$x[i] - u$x_range[i], u$y[i] - u$y_range[i],
                 u$x[i] + u$x_range[i], u$y[i] + u$y_range[i])
}), crs = 3435)
uci_rect_m <- st_transform(uci_rect_ft, 26916)

# CTA-centered SQUARE queries (250/500/750m half-width).
# Squares, not circles, make the UCI-vs-CTA comparison the SAME query shape.
radii <- c(250L, 500L, 750L)
station_idx <- rep(seq_len(nrow(stations)), times = length(radii))
halfwidth_m <- rep(radii, each = nrow(stations))
stations_sf <- st_transform(
  st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326,
           remove = FALSE), 26916)
xy <- st_coordinates(stations_sf)
cta_rect_m <- st_sfc(lapply(seq_along(station_idx), function(i) {
  j <- station_idx[i]
  h <- halfwidth_m[i]
  make_rectangle(xy[j, 1] - h, xy[j, 2] - h,
                 xy[j, 1] + h, xy[j, 2] + h)
}), crs = 26916)

# Reuse 2024 incident records: do NOT use old UCI labels as 2024 truth.
valid_geo <- is.finite(crime$latitude) & is.finite(crime$longitude)
crime <- crime[valid_geo, , drop = FALSE]
crimes_sf <- st_transform(
  st_as_sf(crime, coords = c("longitude", "latitude"), crs = 4326),
  26916)
crime_types <- toupper(trimws(as.character(crime$primary_type)))

queries_m <- c(uci_rect_m, cta_rect_m)
query_source <- c(rep("UCI_original_rectangle", nrow(u)),
                  rep("CTA_centered_rectangle", length(station_idx)))
cat("Performing", length(queries_m), "exact 2024 spatial queries...\n")
incidents_in_query <- st_intersects(queries_m, crimes_sf, sparse = TRUE)
stations_in_query <- st_intersects(queries_m, stations_sf, sparse = TRUE)

count_type <- function(type) {
  vapply(incidents_in_query,
         function(index) as.integer(sum(crime_types[index] == type, na.rm = TRUE)),
         integer(1))
}
station_entries <- as.numeric(stations$total_entries_2024)
cta_entries_in_query <- vapply(stations_in_query,
  function(index) sum(station_entries[index], na.rm = TRUE), numeric(1))

n_uci <- nrow(u)
n_cta <- length(station_idx)
all_counts <- lengths(incidents_in_query)
all_areas_km2 <- as.numeric(st_area(queries_m)) / 1e6
if (any(!is.finite(all_areas_km2)) || any(all_areas_km2 <= 0)) {
  stop("Invalid computed query areas. Inspect source geometry.")
}

comparison <- data.frame(
  query_id = c(paste0("UCI_ROW_", u$source_row_number),
               paste0("CTA_", stations$station_id[station_idx], "_SQ_",
                      halfwidth_m, "M")),
  query_source = query_source,
  geometry_type = "rectangle",
  query_halfwidth_m = c(rep(NA_real_, n_uci), as.numeric(halfwidth_m)),
  station_id = c(rep(NA_character_, n_uci),
                 as.character(stations$station_id[station_idx])),
  station_name = c(rep(NA_character_, n_uci),
                   as.character(stations$station_name[station_idx])),
  area_km2 = all_areas_km2,
  original_uci_count = c(u$count, rep(NA_real_, n_cta)),
  original_uci_sum_arrests = c(u$sum_, rep(NA_real_, n_cta)),
  original_uci_avg_beat = c(u$avg, rep(NA_real_, n_cta)),
  crimes_2024 = all_counts,
  theft_2024 = count_type("THEFT"),
  robbery_2024 = count_type("ROBBERY"),
  battery_2024 = count_type("BATTERY"),
  assault_2024 = count_type("ASSAULT"),
  stations_inside_query = lengths(stations_in_query),
  total_cta_entries_at_stations_inside_query_2024 = cta_entries_in_query,
  center_lon = c(lonlat[, 1], stations$longitude[station_idx]),
  center_lat = c(lonlat[, 2], stations$latitude[station_idx]),
  stringsAsFactors = FALSE
)
comparison$crimes_per_km2_2024 <- comparison$crimes_2024 / comparison$area_km2
# NOTE: CTA entry totals are station ENTRY sums, not human exposure within
# the query's full land area. Rates are intentionally not derived for UCI
# rectangles, which may contain multiple stations or none.

out_comparison <- file.path(processed_dir,
                            "04b_spatial_query_workloads_2024.csv")
out_enriched_uci <- file.path(processed_dir,
                              "04b_uci_range_queries_enriched_2024.csv")
write.csv(comparison, out_comparison, row.names = FALSE)
write.csv(comparison[comparison$query_source == "UCI_original_rectangle", ],
          out_enriched_uci, row.names = FALSE)

# Save sampling metadata so the selected UCI queries can be reproduced.
write.csv(data.frame(source_row_number = u$source_row_number,
                     center_x_stateplane_ft = u$x,
                     center_y_stateplane_ft = u$y,
                     halfwidth_x_ft = u$x_range,
                     halfheight_y_ft = u$y_range,
                     uci_original_count = u$count,
                     query_area_km2 = comparison$area_km2[seq_len(n_uci)]),
          file.path(out_dir, "04b_uci_query_sample.csv"), row.names = FALSE)

# ---- Baseline distributions: preserve query-source labels -------------
plot_group <- ifelse(comparison$query_source == "UCI_original_rectangle",
                     "UCI: sampled rectangles",
                     paste0("CTA: squares +/-",
                            comparison$query_halfwidth_m, " m"))
plot_group <- factor(plot_group, levels = c("UCI: sampled rectangles",
  "CTA: squares +/-250 m", "CTA: squares +/-500 m", "CTA: squares +/-750 m"))

png(file.path(out_dir, "04b_query_area_comparison.png"),
    width = 1350, height = 760, res = 125)
op <- par(mar = c(10, 5, 4, 2))
boxplot(log10(comparison$area_km2) ~ plot_group, las = 2,
        col = "grey80", outline = FALSE,
        xlab = "", ylab = "log10(query area in km^2)",
        main = "Workload comparison: geographic query area")
par(op)
dev.off()

png(file.path(out_dir, "04b_query_count_comparison.png"),
    width = 1350, height = 760, res = 125)
op <- par(mar = c(10, 5, 4, 2))
boxplot(log10(comparison$crimes_2024 + 1) ~ plot_group, las = 2,
        col = "grey80", outline = FALSE,
        xlab = "", ylab = "log10(2024 incident count + 1)",
        main = "Same-year answers: UCI vs CTA-centered rectangles")
par(op)
dev.off()

# Characterize overlap and geographic distribution before any ML claims.
workload_summary <- aggregate(cbind(area_km2, crimes_2024,
                                    stations_inside_query) ~ plot_group,
  data = transform(comparison, plot_group = plot_group),
  FUN = median)
n_with_stations <- aggregate(comparison$stations_inside_query > 0,
                            by = list(query_group = plot_group), FUN = mean)
names(n_with_stations)[2] <- "fraction_with_cta_station"
report <- c(
  "CSP 571 - THE TRANSIT EXPOSURE PARADOX",
  "STEP 4B - ENRICH ORIGINAL UCI RECTANGULAR QUERIES",
  paste("Generated:", Sys.time()),
  paste("Original UCI rows available:", sum(valid)),
  paste("Sampled UCI rectangles:", n_uci),
  paste("New CTA-centered SQUARE queries:", n_cta),
  paste("Total augmented workload rows:", nrow(comparison)),
  paste("2024 Chicago crime points processed:", nrow(crime)),
  paste("CTA stations used:", nrow(stations)),
  paste("UCI State Plane EPSG:3435 center sanity check:",
        sum(plausible), "of", length(plausible), "plausible"),
  "", "MEDIAN VALUES BY QUERY GROUP:",
  capture.output(print(workload_summary, row.names = FALSE)),
  "", "FRACTION OF QUERIES CONTAINING AT LEAST ONE CTA STATION:",
  capture.output(print(n_with_stations, row.names = FALSE)),
  "", "CRITICAL INTERPRETATION:",
  "- UCI's historical count, SUM, AVG are NOT 2024 labels.",
  "- New crimes_2024 was computed for both groups using identical incident data.",
  "- These groups have different areas; raw count comparisons are not controlled comparisons.",
  "- Our CTA-centered geometry is SQUARE here for consistency with UCI rectangles.",
  "- Previous 250/500/750m CTA CIRCLES remain separate analyses in Step 3/4A.",
  "- UCI circular query coordinates are not georeferenced by this script.",
  "- Overlapping queries share incidents and are not independent observations.",
  "- CTA station entries are not complete nearby pedestrian exposure.",
  "- Some CTA stations may be outside Chicago incident-data coverage.",
  "", "OUTPUT FILES:",
  paste(" -", out_comparison),
  paste(" -", out_enriched_uci),
  " - outputs/04b_uci_query_sample.csv",
  " - outputs/04b_query_area_comparison.png",
  " - outputs/04b_query_count_comparison.png"
)
writeLines(report, file.path(out_dir, "04b_uci_augmentation_summary.txt"))
cat("\nSUCCESS - STEP 4B completed.\n")
cat("Enriched", n_uci, "original UCI queries and created", n_cta,
    "CTA-centered rectangular queries.\n")
cat("See outputs/04b_uci_augmentation_summary.txt\n")
