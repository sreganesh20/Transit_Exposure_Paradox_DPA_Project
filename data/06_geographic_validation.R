library(sf)
library(tidyverse)
library(ggplot2)

sf::sf_use_s2(FALSE)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

boundary_file <- if (file.exists("data/raw/Boundaries_-_City_20261121.geojson")) {
  "data/raw/Boundaries_-_City_20261121.geojson"
} else {
  "data/Boundaries_-_City_20261121.geojson"
}

metric_crs <- 26916

city_metric <- st_read(boundary_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(metric_crs)

w  <- read_csv("data/processed/04b_spatial_query_workloads_2024.csv", show_col_types = FALSE)
cr <- read_csv("data/processed/station_crime_ridership_2024.csv", show_col_types = FALSE)

# Exact 426 CTA-centered square queries
cta_q <- w %>%
  filter(query_source == "CTA_centered_rectangle") %>%
  mutate(row_idx = row_number())

centers_sf <- cta_q %>%
  st_as_sf(coords = c("center_lon", "center_lat"), crs = 4326, remove = FALSE) %>%
  st_transform(metric_crs)

xy <- st_coordinates(centers_sf)
hw <- centers_sf$query_halfwidth_m

make_sq <- function(x, y, h) {
  m <- matrix(
    c(x - h, y - h,
      x + h, y - h,
      x + h, y + h,
      x - h, y + h,
      x - h, y - h),
    ncol = 2, byrow = TRUE
  )
  st_polygon(list(m))
}

sq_sf <- centers_sf %>%
  mutate(
    geometry = st_sfc(pmap(list(xy[, 1], xy[, 2], hw), make_sq), crs = metric_crs),
    square_area_m2 = as.numeric(st_area(geometry))
  ) %>%
  st_as_sf()

overlap_tbl <- st_intersection(sq_sf %>% select(row_idx, square_area_m2), city_metric) %>%
  mutate(inside_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  group_by(row_idx) %>%
  summarise(inside_area_m2 = sum(inside_area_m2), .groups = "drop")

# Deduplicate cr_lookup to 1 row per non-missing station_id
cr_lookup <- cr %>%
  filter(!is.na(station_id)) %>%
  group_by(station_id) %>%
  summarise(
    reported_crimes_2024 = first(reported_crimes_2024),
    query_radius_m = first(query_radius_m),
    other_station_buffers_overlapping = first(other_station_buffers_overlapping),
    .groups = "drop"
  )

cov_df <- sq_sf %>%
  st_drop_geometry() %>%
  left_join(overlap_tbl, by = "row_idx") %>%
  left_join(cr_lookup, by = "station_id") %>%
  mutate(
    inside_area_m2 = replace_na(inside_area_m2, 0),
    pct_inside_city = pmin(pmax((inside_area_m2 / square_area_m2) * 100, 0), 100),
    coverage_status = case_when(
      inside_area_m2 <= 0 | pct_inside_city <= 0 ~ "outside Chicago",
      pct_inside_city >= 99.99 ~ "fully covered",
      TRUE ~ "partially covered"
    ),
    zero_crime_w  = (crimes_2024 == 0),
    zero_crime_cr = (reported_crimes_2024 == 0)
  )

stopifnot(nrow(cov_df) == 426)

# 1. Save CSV deliverable
write_csv(cov_df, "data/processed/cta_geographic_coverage.csv")

cat("\n=== Total CTA Queries:", nrow(cov_df), "===\n")
print(cov_df %>% count(coverage_status) %>% mutate(pct = round(n / sum(n) * 100, 2)))

cat("\n=== Zero-Crime Queries (reported_crimes_2024 == 0) by Coverage Status ===\n")
print(cov_df %>% filter(zero_crime_cr) %>% count(coverage_status) %>% mutate(pct = round(n / sum(n) * 100, 2)))

cat("\n=== Zero-Crime Queries (crimes_2024 == 0 in w) by Coverage Status ===\n")
print(cov_df %>% filter(zero_crime_w) %>% count(coverage_status) %>% mutate(pct = round(n / sum(n) * 100, 2)))

# 2. Save visualization deliverable
p <- ggplot() +
  geom_sf(data = city_metric, fill = "grey96", color = "grey45", linewidth = 0.5) +
  geom_sf(
    data = sq_sf %>% left_join(cov_df %>% select(row_idx, coverage_status, zero_crime_cr), by = "row_idx"),
    aes(fill = coverage_status, color = zero_crime_cr),
    alpha = 0.65,
    linewidth = 0.35
  ) +
  scale_fill_manual(values = c(
    "fully covered" = "#2b8cbe",
    "partially covered" = "#fdae61",
    "outside Chicago" = "#d7191c"
  )) +
  scale_color_manual(
    values = c(`FALSE` = NA, `TRUE` = "black"),
    labels = c(`FALSE` = "Reported crimes > 0", `TRUE` = "Zero crimes (n = 54)"),
    name = "Station Crime Status"
  ) +
  labs(
    title = "CTA-Centered Square Query Coverage Within Chicago Limits",
    subtitle = "426 CTA queries evaluated in EPSG:26916; black outlines mark zero-crime stations",
    fill = "Coverage Status",
    caption = "Task A | Geographic Coverage Validation (feat/geographic-validation)"
  ) +
  theme_minimal()

ggsave("figs/geographic_coverage_validation.png", plot = p, width = 9, height = 7, dpi = 300)

cov_n <- cov_df %>% count(coverage_status) %>% deframe()
z_n   <- cov_df %>% filter(zero_crime_cr) %>% count(coverage_status) %>% deframe()

get_cnt <- function(vec, nm) {
  val <- vec[nm]
  if (length(val) == 0 || is.na(val)) 0L else as.integer(val)
}

summary_md <- c(
  "# Task A: Geographic Coverage Validation Summary",
  "",
  paste0("- **Total CTA square queries:** ", nrow(cov_df)),
  paste0("- **Fully covered:** ", get_cnt(cov_n, "fully covered")),
  paste0("- **Partially covered:** ", get_cnt(cov_n, "partially covered")),
  paste0("- **Outside Chicago:** ", get_cnt(cov_n, "outside Chicago")),
  paste0("- **Zero-crime queries (`reported_crimes_2024 == 0`):** ", sum(cov_df$zero_crime_cr, na.rm = TRUE)),
  paste0("  - Outside Chicago: ", get_cnt(z_n, "outside Chicago")),
  paste0("  - Partially covered: ", get_cnt(z_n, "partially covered")),
  paste0("  - Fully covered: ", get_cnt(z_n, "fully covered")),
  "",
  "## Interpretation",
  "- Queries classified as **outside Chicago** or **partially covered** have missing or clipped Chicago crime coverage where query footprints extend beyond city limits.",
  "- Zero-crime queries that are **fully covered** inside Chicago are not explained by boundary truncation and reflect other factors (such as low crime counts within the query window/radius or station-level differences)."
)
writeLines(summary_md, "outputs/task_a_summary.md")