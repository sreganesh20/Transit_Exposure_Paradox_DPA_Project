# CSP 571 — The Transit Exposure Paradox
# STEP 3: 2024 Chicago Crime -> CTA station-centered spatial queries
#
# --- 0. Libraries and paths --------------------------------------------
if (!requireNamespace('sf', quietly = TRUE)) {
  stop("Missing R package 'sf'. Run install.packages('sf') once in RStudio, then rerun.")
}
library(sf)
raw_dir <- file.path('data', 'raw')
processed_dir <- file.path('data', 'processed')
output_dir <- 'outputs'
for (folder in c(raw_dir, processed_dir, output_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
stations_path <- file.path(processed_dir, 'cta_station_ridership_2024.csv')
crimes_path <- file.path(raw_dir, 'chicago_crimes_2024.csv')
for (path in c(stations_path, crimes_path)) {
  if (!file.exists(path)) stop('File not found: ', path,
                               '\nCheck filenames and your RStudio working directory (getwd()).')
}

normalize_names <- function(d) {
  names(d) <- gsub('(^_+|_+$)', '',
                   gsub('[^a-z0-9]+', '_', tolower(names(d))))
  d
}
need <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, ' is missing: ', paste(missing, collapse = ', '),
                            '\nActual columns: ', paste(names(d), collapse = ', '))
}
cat('Reading station and crime records...\n')
stations <- normalize_names(read.csv(stations_path, stringsAsFactors = FALSE))
crimes <- normalize_names(read.csv(crimes_path, stringsAsFactors = FALSE))
need(stations, c('station_id', 'station_name', 'latitude', 'longitude',
                 'total_entries_2024'), 'CTA station file')
need(crimes, c('date', 'latitude', 'longitude', 'primary_type'), 'Chicago crime file')
cat('Downloaded crime CSV rows:', nrow(crimes), '\n')

# Avoid accidentally analyzing the 1000-row portal preview instead of full data.
if (nrow(crimes) < 100000L) {
  stop('Crime CSV has fewer than 100,000 rows. It may be a preview/export limit. ',
       'Download ALL records filtered to Year = 2024, not the first page.')
}

# --- 1. Validate date and geography ------------------------------------
parse_crime_dates <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(as.Date(NA), length(x))
  iso <- grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}', x)
  us <- grepl('^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}', x)
  if (any(iso)) out[iso] <- as.Date(substr(x[iso], 1, 10), format = '%Y-%m-%d')
  if (any(us)) {
    # Keep just M/D/YYYY, stripping any time or AM/PM suffix.
    us_dates <- sub('^([0-9]{1,2}/[0-9]{1,2}/[0-9]{4}).*$', '\\1', x[us])
    out[us] <- as.Date(us_dates, format = '%m/%d/%Y')
  }
  out
}
crimes$incident_date <- parse_crime_dates(crimes$date)
if (all(is.na(crimes$incident_date))) {
  stop('Could not parse any crime dates. First date examples: ',
       paste(head(crimes$date, 5), collapse = ' | '))
}
crimes$latitude <- suppressWarnings(as.numeric(crimes$latitude))
crimes$longitude <- suppressWarnings(as.numeric(crimes$longitude))
valid_date <- !is.na(crimes$incident_date) &
  crimes$incident_date >= as.Date('2024-01-01') &
  crimes$incident_date <= as.Date('2024-12-31')
valid_geo <- is.finite(crimes$latitude) & is.finite(crimes$longitude) &
  crimes$latitude > 40 & crimes$latitude < 43 &
  crimes$longitude > -89 & crimes$longitude < -86
cat('Rows dated in 2024:', sum(valid_date), '\n')
cat('Rows with valid coordinates:', sum(valid_geo), '\n')
cat('2024 rows removed due to missing/invalid coordinates:', sum(valid_date & !valid_geo), '\n')

clean <- crimes[valid_date & valid_geo, , drop = FALSE]
if (!nrow(clean)) stop('No valid 2024 crime incidents left after cleaning.')
if ('id' %in% names(clean)) {
  ids <- as.character(clean$id)
  repeated <- !is.na(ids) & nzchar(ids) & duplicated(ids)
  if (any(repeated)) {
    warning('Removed ', sum(repeated), ' duplicate reported-crime IDs.')
    clean <- clean[!repeated, , drop = FALSE]
  }
}
clean$primary_type <- toupper(trimws(as.character(clean$primary_type)))
write.csv(clean[, c(intersect('id', names(clean)),
                     'incident_date', 'primary_type', 'latitude', 'longitude'), drop = FALSE],
          file.path(processed_dir, 'chicago_crimes_2024_clean.csv'),
          row.names = FALSE)

# Note: Missing coordinates and records outside Chicago limit geographic coverage.
stations$latitude <- suppressWarnings(as.numeric(stations$latitude))
stations$longitude <- suppressWarnings(as.numeric(stations$longitude))
stations$total_entries_2024 <- suppressWarnings(as.numeric(stations$total_entries_2024))
station_ok <- is.finite(stations$latitude) & is.finite(stations$longitude) &
  stations$latitude > 40 & stations$latitude < 43 &
  stations$longitude > -89 & stations$longitude < -86 &
  is.finite(stations$total_entries_2024) & stations$total_entries_2024 > 0
stations_excluded <- stations[!station_ok, , drop = FALSE]
stations <- stations[station_ok, , drop = FALSE]
if (!nrow(stations)) stop('No valid CTA stations remain after geography/ridership checks.')
if (anyDuplicated(as.character(stations$station_id))) stop('CTA station IDs are not unique.')

# --- 2. Convert coordinates to meters for spatial queries --------------
# EPSG:26916 = NAD83 / UTM zone 16N. Distances use meters.
# This is separate from any original UCI State Plane coordinates.
crimes_sf <- sf::st_transform(
  sf::st_as_sf(clean, coords = c('longitude', 'latitude'), crs = 4326),
  crs = 26916
)
stations_sf <- sf::st_transform(
  sf::st_as_sf(stations, coords = c('longitude', 'latitude'), crs = 4326),
  crs = 26916
)

# --- 3. Count incidents in three station-centered circles --------------
radii_m <- c(250L, 500L, 750L)
result_list <- vector('list', length(radii_m))
for (i in seq_along(radii_m)) {
  radius <- radii_m[i]
  cat('Spatial queries: radius ', radius, ' m ...\n', sep = '')
  circles <- sf::st_buffer(stations_sf, dist = radius)
  matches <- sf::st_intersects(circles, crimes_sf, sparse = TRUE)
  # If a crime is inside two station circles, it is counted in both
  # station-centered queries (as expected). Never sum these query counts
  # to claim a unique citywide incident count.
  total_count <- lengths(matches)
  type_count <- function(label) {
    vapply(matches, function(idx) sum(clean$primary_type[idx] == label, na.rm = TRUE),
           integer(1))
  }
  r <- stations[, c('station_id', 'station_name', 'latitude', 'longitude',
                    'total_entries_2024', intersect('days_reported', names(stations))),
                drop = FALSE]
  r$query_radius_m <- radius
  r$reported_crimes_2024 <- total_count
  r$theft_2024 <- type_count('THEFT')
  r$robbery_2024 <- type_count('ROBBERY')
  r$battery_2024 <- type_count('BATTERY')
  r$assault_2024 <- type_count('ASSAULT')
  r$crimes_per_100k_entries <- 1e5 * r$reported_crimes_2024 / r$total_entries_2024
  r$raw_count_rank <- rank(-r$reported_crimes_2024, ties.method = 'min')
  r$activity_adjusted_rank <- rank(-r$crimes_per_100k_entries, ties.method = 'min')
  r$rank_change <- r$raw_count_rank - r$activity_adjusted_rank
  # An additional station has overlapping geography if its center is within 2*radius.
  station_overlaps <- sf::st_intersects(circles, sparse = TRUE)
  r$other_station_buffers_overlapping <- pmax(0L, lengths(station_overlaps) - 1L)
  result_list[[i]] <- r
}
final <- do.call(rbind, result_list)
rownames(final) <- NULL
write.csv(final,
          file.path(processed_dir, 'station_crime_ridership_2024.csv'),
          row.names = FALSE)

# --- 4. Inspect 500m rankings and create a chart ------------------------
r500 <- final[final$query_radius_m == 500, , drop = FALSE]
r500 <- r500[order(r500$raw_count_rank, -r500$total_entries_2024), ]
write.csv(r500,
          file.path(output_dir, 'station_rank_shift_500m_2024.csv'),
          row.names = FALSE)

top <- head(r500, 15)
if (nrow(top) > 0) {
  nm <- as.character(top$station_name)
  nm[is.na(nm) | !nzchar(nm)] <- as.character(top$station_id[is.na(nm) | !nzchar(nm)])
  png(file.path(output_dir, 'station_crime_500m_2024.png'),
      width = 1150, height = 820, res = 120)
  op <- par(mar = c(5, 12, 4, 2))
  barplot(rev(top$reported_crimes_2024),
          names.arg = rev(nm), horiz = TRUE, las = 1,
          col = 'grey55', border = NA,
          main = '2024 reported crimes within 500 m of CTA stations',
          xlab = 'Reported incidents in 500 m query (NOT transit-only crimes)')
  par(op)
  dev.off()
}

# --- Additional plots: 250m and 750m ----------------------------

for (radius in c(250, 750)) {
  
  # Filter results for the selected radius
  r_data <- final[final$query_radius_m == radius, ]
  
  # Sort stations by reported crime count
  r_data <- r_data[
    order(-r_data$reported_crimes_2024),
  ]
  
  # Distinguish stations sharing the same name
  labels <- as.character(r_data$station_name)
  duplicate_names <- duplicated(labels) |
    duplicated(labels, fromLast = TRUE)
  
  labels[duplicate_names] <- paste0(
    labels[duplicate_names],
    " (", r_data$station_id[duplicate_names], ")"
  )
  
  # Select top 15 stations
  top <- head(r_data, 15)
  top_labels <- head(labels, 15)
  
  # Save plot
  png(
    file.path(
      output_dir,
      paste0("station_crime_", radius, "m_2024.png")
    ),
    width = 1300,
    height = 850,
    res = 120
  )
  
  op <- par(mar = c(5, 15, 4, 2))
  
  barplot(
    rev(top$reported_crimes_2024),
    names.arg = rev(top_labels),
    horiz = TRUE,
    las = 1,
    col = "grey55",
    border = NA,
    main = paste0(
      "2024 Reported Crimes Within ",
      radius,
      "m of CTA Stations"
    ),
    xlab = "Reported incidents (NOT transit-only crimes)"
  )
  
  par(op)
  dev.off()
  
  cat("Generated plot for", radius, "meters\n")
}

# Rank correlations describe whether station orders differ by metric.
rho <- suppressWarnings(cor(r500$raw_count_rank,
                             r500$activity_adjusted_rank,
                             method = 'spearman', use = 'complete.obs'))
summary_lines <- c(
  'CSP 571 - THE TRANSIT EXPOSURE PARADOX',
  'STEP 3: 2024 CRIME + CTA STATION SPATIAL QUERY JOIN',
  paste('Generated:', Sys.time()),
  paste('Original downloaded crime records:', nrow(crimes)),
  paste('2024 crimes before coordinate validation:', sum(valid_date)),
  paste('2024 crimes with valid coordinates after ID deduplication:', nrow(clean)),
  paste('2024 crimes with missing/invalid geography:', sum(valid_date & !valid_geo)),
  paste('CTA stations processed:', nrow(stations)),
  paste('CTA stations excluded due to invalid geography/entries:', nrow(stations_excluded)),
  paste('Radii (meters):', paste(radii_m, collapse = ', ')),
  paste('Spatial query records generated:', nrow(final)),
  paste('Stations with overlapping 500m buffers:',
        sum(r500$other_station_buffers_overlapping > 0)),
  paste('Spearman correlation, 500m raw-count vs activity-adjusted ranks:',
        signif(rho, 4)),
  '',
  'TOP 15 STATIONS BY 500M RAW REPORTED-CRIME COUNT:',
  capture.output(print(head(r500[, c('station_name', 'reported_crimes_2024',
                                     'total_entries_2024', 'crimes_per_100k_entries',
                                     'raw_count_rank', 'activity_adjusted_rank')], 15),
                       row.names = FALSE)),
  '',
  'INTERPRETATION WARNINGS:',
  '- Nearby crime records are NOT proof that a crime occurred on CTA property.',
  '- CTA station entries are NOT total pedestrian exposure or individual risk.',
  '- Crime locations are block-level approximate; 250m circles are sensitive to this.',
  '- Overlapping buffers intentionally share incidents. Never sum station counts as unique crimes.',
  '- Station coordinates may reflect a newer station network than in 2024.',
  '- Stations beyond Chicago city limits and boundary-truncated buffers are not filtered yet.',
  '- Original UCI query geometries will be integrated in a later stage.'
)
writeLines(summary_lines, file.path(output_dir, 'crime_spatial_summary_2024.txt'))
cat('\nSUCCESS - Step 3 spatial crime join completed.\n')
cat('2024 crime records with valid locations:', nrow(clean), '\n')
cat('Stations:', nrow(stations), '\n')
cat('Station-centered queries:', nrow(final), '\n')
cat('Saved:\n')
cat('  data/processed/station_crime_ridership_2024.csv\n')
cat('  data/processed/chicago_crimes_2024_clean.csv\n')
cat('  outputs/crime_spatial_summary_2024.txt\n')
cat('  outputs/station_rank_shift_500m_2024.csv\n')
cat('  outputs/station_crime_500m_2024.png\n')