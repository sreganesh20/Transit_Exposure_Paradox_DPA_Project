# CSP 571 – The Transit Exposure Paradox
# Step 5B: Check transfer evaluation for zero-count/coverage artifacts.
# Uses existing Step 2 and Step 5 files. No new downloads or packages.
# NOTE: Filtering zero-count queries is a DIAGNOSTIC sensitivity analysis,
# not a valid final city-boundary filter. Do not treat zero counts as errors
# without confirming each query's geographic coverage.

pred_path <- file.path("data", "processed", "05_model_predictions_2024.csv")
station_path <- file.path("data", "processed", "cta_station_ridership_2024.csv")
out_dir <- "outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(pred_path, station_path)) {
  if (!file.exists(f)) stop("Missing ", f, ". Run previous steps first.")
}
x <- read.csv(pred_path, stringsAsFactors = FALSE)
stations <- read.csv(station_path, stringsAsFactors = FALSE)
needed <- c("query_source", "station_id", "query_halfwidth_m",
            "crimes_2024", "abs_log_error_baseline",
            "abs_log_error_augmented", "within_uci_training_support")
if (!all(needed %in% names(x))) {
  stop("Step 5 file missing: ", paste(setdiff(needed, names(x)), collapse = ", "))
}
if (!all(c("station_id", "station_name", "latitude", "longitude") %in% names(stations))) {
  stop("Unexpected CTA station file columns.")
}

cta <- subset(x, query_source == "CTA_centered_rectangle")
uci <- subset(x, evaluation_set == "UCI_spatial_holdout")
if (nrow(cta) != 426L) {
  warning("Expected 426 CTA queries from Step 4B, found ", nrow(cta))
}
if (anyNA(cta$station_id)) stop("CTA station IDs missing.")
cta$station_id <- as.character(cta$station_id)
stations$station_id <- as.character(stations$station_id)
cta$within_training_support <- as.character(cta$within_uci_training_support) == "TRUE"
cta$zero_incident_count <- cta$crimes_2024 == 0
cta$log_error_gain <- cta$abs_log_error_baseline - cta$abs_log_error_augmented
cta$radius <- paste0(as.integer(cta$query_halfwidth_m), " m halfwidth")

# Each item below is a descriptive subgroup, not an independent test set.
subgroups <- list(
  CTA_all = rep(TRUE, nrow(cta)),
  CTA_within_training_support = cta$within_training_support,
  CTA_outside_training_support = !cta$within_training_support,
  CTA_zero_recorded_incidents = cta$zero_incident_count,
  CTA_nonzero_recorded_incidents = !cta$zero_incident_count,
  CTA_supported_nonzero = cta$within_training_support & !cta$zero_incident_count
)
for (r in c(250, 500, 750)) {
  subgroups[[paste0("CTA_", r, "m_square")]] <- cta$query_halfwidth_m == r
}

score <- function(df, group) {
  n <- nrow(df)
  if (!n) return(data.frame(
    group = group, queries = 0L, stations = 0L,
    mean_crimes = NA_real_, baseline_mae_log = NA_real_,
    augmented_mae_log = NA_real_, baseline_minus_augmented = NA_real_,
    reduction_pct = NA_real_
  ))
  b <- mean(df$abs_log_error_baseline)
  a <- mean(df$abs_log_error_augmented)
  data.frame(
    group = group,
    queries = n,
    stations = length(unique(df$station_id)),
    mean_crimes = mean(df$crimes_2024),
    baseline_mae_log = b,
    augmented_mae_log = a,
    baseline_minus_augmented = b - a,
    reduction_pct = 100 * (b - a) / b
  )
}
breakdown <- do.call(rbind, Map(
  function(name, flag) score(cta[flag, , drop = FALSE], name),
  names(subgroups), subgroups
))
rownames(breakdown) <- NULL
write.csv(breakdown, file.path(out_dir, "05b_transfer_error_breakdown.csv"),
          row.names = FALSE)

# Audit geographically suspicious queries: 0 records in a whole
# station-centered query may indicate an out-of-jurisdiction station.
# This CSV is for verification against the official Chicago city polygon.
zero_q <- cta[cta$zero_incident_count, c("station_id", "query_halfwidth_m",
                  "center_lat", "center_lon", "crimes_2024",
                  "abs_log_error_baseline", "abs_log_error_augmented",
                  "within_training_support"), drop = FALSE]
station_lookup <- unique(stations[, c("station_id", "station_name",
                                     "latitude", "longitude")])
if (anyDuplicated(station_lookup$station_id)) {
  stop("Station geographic lookup contains duplicate IDs.")
}
zero_q <- merge(zero_q, station_lookup, by = "station_id", all.x = TRUE)
zero_q <- zero_q[order(zero_q$station_name, zero_q$query_halfwidth_m), ]
write.csv(zero_q, file.path(out_dir, "05b_zero_count_queries_to_verify.csv"),
          row.names = FALSE)

# Contrast the same metric on the untouched UCI holdout.
uci_b <- mean(uci$abs_log_error_baseline, na.rm = TRUE)
uci_a <- mean(uci$abs_log_error_augmented, na.rm = TRUE)

# Plot: zero vs nonzero and inside vs outside training support.
plot_groups <- c("CTA_within_training_support", "CTA_outside_training_support",
                 "CTA_zero_recorded_incidents", "CTA_nonzero_recorded_incidents")
z <- breakdown[match(plot_groups, breakdown$group), ]
mat <- rbind(z$baseline_mae_log, z$augmented_mae_log)
colnames(mat) <- c("In UCI\nfeature range", "Outside UCI\nfeature range",
                   "Zero recorded\ncrimes", "Nonzero\ncrimes")
rownames(mat) <- c("Baseline", "Augmented")
png(file.path(out_dir, "05b_transfer_diagnostics.png"),
    width = 1250, height = 850, res = 125)
op <- par(mar = c(7, 5, 4, 2))
barplot(mat, beside = TRUE, col = c("grey40", "grey75"),
        las = 1, ylim = c(0, max(mat, na.rm = TRUE) * 1.16),
        ylab = "Mean absolute log(1 + crime count) error",
        main = "Transfer diagnostics: where does CTA query error come from?",
        legend.text = rownames(mat),
        args.legend = list(x = "topright", bty = "n"))
par(op)
dev.off()

summary_text <- c(
  "CSP 571 – TRANSIT EXPOSURE PARADOX",
  "STEP 5B: COVERAGE / EXTRAPOLATION DIAGNOSTICS",
  paste("Generated:", Sys.time()),
  "",
  paste("UCI holdout baseline MAE (log):", signif(uci_b, 5)),
  paste("UCI holdout augmented MAE (log):", signif(uci_a, 5)),
  "",
  "CTA SUBGROUP METRICS:",
  capture.output(print(breakdown, row.names = FALSE, digits = 4)),
  "",
  paste("Distinct station IDs with >=1 zero-count query:",
        length(unique(zero_q$station_id))),
  paste("Zero-count queries:", nrow(zero_q), "of", nrow(cta)),
  paste("CTA queries in UCI central feature range:",
        sum(cta$within_training_support), "of", nrow(cta)),
  "",
  "WARNINGS:",
  "- CTA stations outside Chicago city limits are NOT covered by Chicago Police Department data.",
  "- A query with 0 crimes is NOT evidence a station is crime-free.",
  "- Some Chicago city-edge buffers may contain land outside city jurisdiction.",
  "- Excluding zero counts is ONLY a sensitivity check; official boundary clipping is the next required validation.",
  "- Within-support uses marginal central ranges for area, longitude and latitude, not a causal or true multivariate support test.",
  "- A large transfer error may reflect query geography/size, jurisdiction, missing records, and model misspecification.",
  "- Do not interpret ridership as passenger-level crime exposure."
)
writeLines(summary_text, file.path(out_dir, "05b_transfer_diagnostics_summary.txt"))
cat("\nSUCCESS - Step 5B diagnostic complete.\n")
cat("Saved diagnostics, zero-count query list, plot, and summary in outputs/.\n")
print(breakdown, row.names = FALSE)
