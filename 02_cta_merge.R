# CSP 571 — The Transit Exposure Paradox
# Step 2: Join CTA daily station ridership (2024) to station geography.
#
# ---- 1. Locations and input checks -----------------------------------
raw_dir <- file.path('data', 'raw')
processed_dir <- file.path('data', 'processed')
output_dir <- 'outputs'
for (folder in c(raw_dir, processed_dir, output_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
ridership_path <- file.path(raw_dir, 'cta_ridership.csv')
stops_path <- file.path(raw_dir, 'cta_stops.csv')
for (path in c(ridership_path, stops_path)) {
  if (!file.exists(path)) {
    stop('MISSING FILE: ', path,
         '\nDownload the official CTA CSVs and use the filenames above.')
  }
}

cat('Reading CTA ridership and station-location CSV files...\n')
rides <- read.csv(ridership_path, stringsAsFactors = FALSE,
                  check.names = FALSE)
stops <- read.csv(stops_path, stringsAsFactors = FALSE,
                  check.names = FALSE)

normalize_names <- function(d) {
  names(d) <- gsub('(^_+|_+$)', '',
                   gsub('[^a-z0-9]+', '_', tolower(names(d))))
  d
}
rides <- normalize_names(rides)
stops <- normalize_names(stops)
cat('Ridership columns: ', paste(names(rides), collapse = ', '), '\n')
cat('Stops columns: ', paste(names(stops), collapse = ', '), '\n')

require_columns <- function(d, fields, what) {
  missing <- setdiff(fields, names(d))
  if (length(missing) > 0L) {
    stop(what, ' is missing columns: ', paste(missing, collapse = ', '))
  }
}
require_columns(rides, c('station_id', 'stationname', 'date', 'rides'),
                'Ridership file')
require_columns(stops, c('map_id', 'station_name'), 'Stops file')

# ---- 2. Parse daily records and isolate calendar year 2024 ------------
parse_cta_dates <- function(x) {
  x <- trimws(as.character(x))
  # Both Socrata ISO timestamps and US-formatted dates are supported.
  datepart <- sub('[ T].*$', '', x)
  result <- as.Date(datepart, format = '%Y-%m-%d')
  needs_us <- is.na(result) & nzchar(datepart)
  result[needs_us] <- as.Date(datepart[needs_us], format = '%m/%d/%Y')
  result
}
rides$date <- parse_cta_dates(rides$date)
if (anyNA(rides$date)) {
  stop('Some ridership dates could not be parsed. Inspect the date column.')
}
rides$station_id <- trimws(as.character(rides$station_id))
rides$rides <- suppressWarnings(as.numeric(gsub(',', '',
                                                as.character(rides$rides))))
rides_2024 <- rides[rides$date >= as.Date('2024-01-01') &
                    rides$date <= as.Date('2024-12-31'),
                    c('station_id', 'stationname', 'date', 'rides'), drop = FALSE]
if (nrow(rides_2024) == 0L) {
  stop('No 2024 records found. Check whether the downloaded ridership file includes 2024.')
}
if (anyNA(rides_2024$rides) || any(rides_2024$rides < 0)) {
  stop('Missing or negative 2024 ridership values; examine the source before proceeding.')
}

# There should be one daily record per station. Never silently double count.
key <- paste(rides_2024$station_id, rides_2024$date, sep = '|')
if (anyDuplicated(key)) {
  write.csv(rides_2024[duplicated(key) | duplicated(key, fromLast = TRUE), ],
            file.path(output_dir, 'cta_duplicate_station_dates.csv'),
            row.names = FALSE)
  stop('Duplicate station/date records found. See outputs/cta_duplicate_station_dates.csv')
}
cat('2024 daily station-entry records:', format(nrow(rides_2024), big.mark = ','), '\n')

# ---- 3. Extract stop coordinates --------------------------------------
# The official stops file typically has Location = '(latitude, longitude)'.
# Alternative exports may have separate latitude/longitude columns.
extract_point <- function(s) {
  tokens <- regmatches(as.character(s),
                       gregexpr('-?[0-9]+(?:\\.[0-9]+)?', as.character(s),
                                perl = TRUE))[[1]]
  nums <- suppressWarnings(as.numeric(tokens))
  lat <- nums[nums >= 40 & nums <= 43]
  lon <- nums[nums >= -89 & nums <= -86]
  if (length(lat) < 1L || length(lon) < 1L) return(c(NA_real_, NA_real_))
  c(lat[1], lon[1])
}

if (all(c('latitude', 'longitude') %in% names(stops))) {
  stops$latitude <- suppressWarnings(as.numeric(stops$latitude))
  stops$longitude <- suppressWarnings(as.numeric(stops$longitude))
} else if ('location' %in% names(stops)) {
  point_list <- lapply(stops$location, extract_point)
  point_matrix <- do.call(rbind, point_list)
  stops$latitude <- point_matrix[, 1]
  stops$longitude <- point_matrix[, 2]
} else {
  stop('No coordinate fields found. Expected Location or latitude/longitude.')
}
stops$station_id <- trimws(as.character(stops$map_id))
stops$station_name <- trimws(as.character(stops$station_name))
stops <- stops[nzchar(stops$station_id) & !is.na(stops$station_id), ]

# Stop file contains multiple platforms/stops per parent station (MAP_ID).
# Use a median position for each parent station, rather than duplicating rides.
valid <- is.finite(stops$latitude) & is.finite(stops$longitude) &
  stops$latitude >= 40 & stops$latitude <= 43 &
  stops$longitude >= -89 & stops$longitude <= -86
stops_valid <- stops[valid, ]
if (nrow(stops_valid) == 0L) {
  stop('No valid geographic coordinates parsed from the CTA stops file.')
}
station_geography <- aggregate(
  cbind(latitude, longitude) ~ station_id,
  data = stops_valid, FUN = median
)
station_geography$station_name <- vapply(
  station_geography$station_id,
  function(id) {
    names_at_id <- stops_valid$station_name[stops_valid$station_id == id]
    names_at_id <- names_at_id[!is.na(names_at_id) & nzchar(names_at_id)]
    if (length(names_at_id)) names_at_id[1] else NA_character_
  }, character(1)
)

# ---- 4. Annual summaries + ID-based geographic join -------------------
annual <- aggregate(rides ~ station_id, data = rides_2024, FUN = sum)
names(annual)[names(annual) == 'rides'] <- 'total_entries_2024'
day_counts <- aggregate(date ~ station_id, data = rides_2024, FUN = length)
names(day_counts)[names(day_counts) == 'date'] <- 'days_reported'
annual <- merge(annual, day_counts, by = 'station_id', all.x = TRUE)

# Join via station_id = MAP_ID (NOT by potentially ambiguous station names).
joined <- merge(annual, station_geography, by = 'station_id', all.x = TRUE)
joined <- joined[order(-joined$total_entries_2024), ]
rownames(joined) <- NULL

matched <- is.finite(joined$latitude) & is.finite(joined$longitude)
missing_coordinates <- joined[!matched, ]
write.csv(missing_coordinates,
          file.path(output_dir, 'cta_unmatched_stations_2024.csv'),
          row.names = FALSE)
write.csv(joined,
          file.path(processed_dir, 'cta_station_ridership_2024.csv'),
          row.names = FALSE)
write.csv(rides_2024,
          file.path(processed_dir, 'cta_daily_ridership_2024.csv'),
          row.names = FALSE)

# ---- 5. Baseline output and chart -------------------------------------
summary_lines <- c(
  'CSP 571 - THE TRANSIT EXPOSURE PARADOX',
  'STEP 2: CTA RIDERSHIP + STATION GEOGRAPHY',
  paste('Generated:', Sys.time()),
  'Period: calendar year 2024 (366 days)',
  paste('Daily station-entry rows:', nrow(rides_2024)),
  paste('Unique station IDs in 2024 ridership:', nrow(annual)),
  paste('Parent station IDs with valid location:', nrow(station_geography)),
  paste('Ridership stations matched to geography:', sum(matched)),
  paste('Ridership stations WITHOUT valid geography:', sum(!matched)),
  paste('Stations with fewer than 366 reporting days:',
        sum(joined$days_reported < 366)),
  paste('Total entries across all ridership stations:',
        format(sum(joined$total_entries_2024), big.mark = ',')),
  'NOTE: Entries count turnstile entries, not all people in nearby streets.',
  'NOTE: Station location file may reflect a newer system than 2024.',
  'NOTE: Some CTA stations are outside Chicago; geographic filtering comes later.',
  '',
  'TOP 15 STATIONS BY 2024 ENTRIES:',
  capture.output(print(head(joined[c('station_id', 'station_name',
                                    'total_entries_2024', 'days_reported')], 15),
                       row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, 'cta_summary_2024.txt'))

# Plot highest-entry matched stations, avoiding NA labels.
top <- head(joined[matched, ], 15)
if (nrow(top) > 0L) {
  labels <- ifelse(is.na(top$station_name) | !nzchar(top$station_name),
                   top$station_id, top$station_name)
  png(file.path(output_dir, 'cta_top15_ridership_2024.png'),
      width = 1150, height = 780, res = 120)
  op <- par(mar = c(5, 12, 4, 2))
  barplot(rev(top$total_entries_2024 / 1e6),
          names.arg = rev(labels), horiz = TRUE, las = 1,
          xlab = 'Total 2024 CTA station entries (millions)',
          main = 'Highest-entry CTA rail stations, 2024',
          col = 'grey60', border = NA)
  par(op)
  dev.off()
}

cat('\nSUCCESS — Step 2 CTA join completed.\n')
cat('Stations in 2024 ridership:', nrow(joined), '\n')
cat('Matched to locations:', sum(matched), '\n')
cat('Unmatched:', sum(!matched), '\n')
cat('Saved:\n')
cat('  data/processed/cta_station_ridership_2024.csv\n')
cat('  data/processed/cta_daily_ridership_2024.csv\n')
cat('  outputs/cta_summary_2024.txt\n')
cat('  outputs/cta_unmatched_stations_2024.csv\n')
cat('  outputs/cta_top15_ridership_2024.png\n')
