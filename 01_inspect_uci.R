# CSP 571 — The Transit Exposure Paradox

# ---- 1. Create local project folders ---------------------------------
folders <- c('data/raw', 'data/processed', 'outputs', 'scripts')
for (folder in folders) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Get the original UCI ZIP --------------------------------------
uci_url <- paste0(
  'https://archive.ics.uci.edu/static/public/493/',
  'query%2Banalytics%2Bworkloads%2Bdataset.zip'
)
zip_path <- file.path('data', 'raw', 'uci_query_workloads.zip')

if (!file.exists(zip_path) || file.size(zip_path) < 10000) {
  cat('Downloading the official UCI dataset...\n')
  tryCatch(
    utils::download.file(uci_url, destfile = zip_path, mode = 'wb', quiet = FALSE),
    error = function(e) {
      stop('UCI download failed. Download the ZIP manually from ',
           'https://archive.ics.uci.edu/dataset/493/query%2Banalytics%2Bworkloads%2Bdataset',
           ' and save it as data/raw/uci_query_workloads.zip. Details: ', conditionMessage(e))
    }
  )
}

out_dir <- file.path('data', 'raw', 'uci')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::unzip(zip_path, exdir = out_dir)

# ---- 3. Read all three published CSVs --------------------------------
all_csvs <- list.files(out_dir, pattern = '\\.csv$', full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)
expected_files <- c('Radius-Queries.csv', 'Radius-Queries-Count.csv',
                    'Range-Queries-Aggregates.csv')

find_csv <- function(filename) {
  hits <- all_csvs[tolower(basename(all_csvs)) == tolower(filename)]
  if (length(hits) != 1L) stop('Expected exactly one ', filename,
                              '. Found ', length(hits), '. Check the ZIP contents.')
  hits
}

expected_columns <- list(
  'Radius-Queries.csv' = c('X-coordinate', 'Y-coordinate', 'R-Radius'),
  'Radius-Queries-Count.csv' = c('X-coordinate', 'Y-coordinate', 'R-Radius', 'Count'),
  'Range-Queries-Aggregates.csv' = c('X-coordinate', 'Y-coordinate',
                                      'X-range', 'Y-range', 'Count', 'SUM', 'AVG')
)

read_uci_csv <- function(path, expected) {
  line1 <- readLines(path, n = 1, warn = FALSE)
  parts <- trimws(gsub('"', '', strsplit(line1, ',', fixed = TRUE)[[1]]))
  numeric_first_line <- length(parts) > 0L &&
    all(!is.na(suppressWarnings(as.numeric(parts))))
  d <- utils::read.csv(path, header = !numeric_first_line,
                       check.names = FALSE, stringsAsFactors = FALSE)
  # Some UCI files have no header; assign names based on official documentation.
  if (numeric_first_line) {
    if (ncol(d) == length(expected)) {
      names(d) <- expected
    } else if (ncol(d) == length(expected) + 1L) {
      names(d) <- c('ID', expected)
    } else {
      warning('Unexpected column count in ', basename(path), ': ', ncol(d))
    }
  }
  d
}

datasets <- setNames(vector('list', length(expected_files)), expected_files)
for (nm in expected_files) {
  datasets[[nm]] <- read_uci_csv(find_csv(nm), expected_columns[[nm]])
}

# ---- 4. Summarize and save reproducible baseline outputs ---------------
get_count <- function(d) {
  normalized <- gsub('[^a-z0-9]', '', tolower(names(d)))
  idx <- which(normalized == 'count')
  if (length(idx) != 1L) return(NULL)
  suppressWarnings(as.numeric(d[[idx]]))
}

inventory <- data.frame(
  file = names(datasets),
  n_rows = vapply(datasets, nrow, integer(1)),
  n_columns = vapply(datasets, ncol, integer(1)),
  columns = vapply(datasets, function(d) paste(names(d), collapse = ' | '), character(1)),
  stringsAsFactors = FALSE
)
write.csv(inventory, 'outputs/uci_inventory.csv', row.names = FALSE)

report <- c(
  'CSP 571 - THE TRANSIT EXPOSURE PARADOX',
  'STEP 1: UCI DATASET BASELINE INVENTORY',
  paste('Generated:', Sys.time()),
  'Source: https://archive.ics.uci.edu/dataset/493/query%2Banalytics%2Bworkloads%2Bdataset',
  'NOTE: Rows are spatial queries, NOT individual crime incidents.',
  ''
)
for (nm in names(datasets)) {
  d <- datasets[[nm]]
  count <- get_count(d)
  report <- c(report, paste0('=== ', nm, ' ==='),
              paste('Rows:', nrow(d), '| Columns:', ncol(d)),
              paste('Column names:', paste(names(d), collapse = ', ')),
              'First three rows:', capture.output(print(utils::head(d, 3))),
              'Summary:', capture.output(summary(d)), '')
  if (!is.null(count)) {
    report <- c(report,
                paste('Count zero fraction:', round(mean(count == 0, na.rm = TRUE), 4)),
                paste('Count median:', stats::median(count, na.rm = TRUE)), '')
  }
}
writeLines(report, 'outputs/uci_overview.txt')

# Plot the distribution of query counts (where present).
count_datasets <- datasets[vapply(datasets, function(d) !is.null(get_count(d)), logical(1))]
if (length(count_datasets) > 0L) {
  png('outputs/uci_query_count_distributions.png', width = 1200, height = 650, res = 130)
  oldpar <- par(mfrow = c(1, length(count_datasets)), mar = c(5, 4, 3, 1))
  for (nm in names(count_datasets)) {
    vals <- get_count(count_datasets[[nm]])
    hist(log10(vals[is.finite(vals) & vals >= 0] + 1), breaks = 40,
         main = nm, xlab = 'log10(Count + 1)', col = 'grey70', border = 'white')
  }
  par(oldpar)
  dev.off()
}

cat('\nSUCCESS — Step 1 inventory completed.\n')
print(inventory, row.names = FALSE)
cat('\nSaved:\n')
cat('  outputs/uci_inventory.csv\n')
cat('  outputs/uci_overview.txt\n')
cat('  outputs/uci_query_count_distributions.png\n')