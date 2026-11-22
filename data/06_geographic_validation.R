# ==============================================================================
# Task A: Geographic Coverage Validation for 426 CTA-Centered Square Queries
#
# Required Packages:
#   sf, tidyverse, ggplot2
#
# Boundary Data Source:
#   Official City of Chicago City Boundaries GeoJSON
#   (e.g., data/raw/Boundaries_-_City_20261121.geojson or
#    data/Boundaries_-_City_20261121.geojson)
# ==============================================================================

library(sf)
library(tidyverse)
library(ggplot2)

sf::sf_use_s2(FALSE)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

boundary_candidates <- c(
  "data/raw/Boundaries_-_City_20261121.geojson",
  "data/Boundaries_-_City_20261121.geojson"
)
boundary_file <- boundary_candidates[file.exists(boundary_candidates)][1]

if (is.na(boundary_file) || !nzchar(boundary_file)) {
  stop(
    "Chicago city boundary GeoJSON not found. Expected at ",
    paste(boundary_candidates, collapse = " or ")
  )
}

metric_crs <- 26916

city_metric <- st_read(boundary_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(metric_crs)

w <- read_csv("data/processed/04b_spatial_query_workloads_2024.csv", show_col_types = FALSE)

# Exact 426 CTA-centered square queries
cta_q <- w %>%
  filter(query_source == "CTA_centered_rectangle") %>%
  mutate(row_idx = row_number())

stopifnot(nrow(cta_q) == 426)
if ("query_id" %in% names(cta_q)) {
  stopifnot(n_distinct(cta_q$query_id) == nrow(cta_q))
}

centers_sf <- cta_q %>%
  st_as_sf(coords = c("center_lon", "center_lat"), crs = 4326, remove = FALSE) %>%
  st_transform(metric_crs)

xy <- st_coordinates(centers_sf)
hw <- centers_sf$query_halfwidth_m

make_sq <- function(x, y, h) {
  m <- matrix(
    c(
      x - h, y - h,
      x + h, y - h,
      x + h, y + h,
      x - h, y + h,
      x - h, y - h
    ),
    ncol = 2,
    byrow = TRUE
  )
  st_polygon(list(m))
}

sq_sf <- centers_sf %>%
  mutate(
    geometry = st_sfc(pmap(list(xy[, 1], xy[, 2], hw), make_sq), crs = metric_crs),
    square_area_m2 = as.numeric(st_area(geometry))
  )

# Intersect square footprints with Chicago boundary in projected CRS (EPSG:26916)
overlap_sf <- suppressWarnings(
  st_intersection(sq_sf %>% select(row_idx, square_area_m2), city_metric)
)

overlap_tbl <- overlap_sf %>%
  mutate(inside_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  group_by(row_idx) %>%
  summarise(
    inside_area_m2 = sum(inside_area_m2, na.rm = TRUE),
    .groups = "drop"
  )

cov_df <- sq_sf %>%
  st_drop_geometry() %>%
  left_join(overlap_tbl, by = "row_idx") %>%
  mutate(
    inside_area_m2 = coalesce(inside_area_m2, 0),
    pct_inside_chicago = pmin(pmax((inside_area_m2 / square_area_m2) * 100, 0), 100),
    coverage_status = case_when(
      pct_inside_chicago >= 99.99 ~ "fully covered",
      pct_inside_chicago <= 0.0001 ~ "outside Chicago",
      TRUE ~ "partially covered"
    ),
    zero_crime_query = (crimes_2024 == 0)
  )

stopifnot(nrow(cov_df) == 426)
if ("query_id" %in% names(cov_df)) {
  stopifnot(n_distinct(cov_df$query_id) == nrow(cov_df))
}

# 1. Save CSV deliverable
write_csv(cov_df, "data/processed/cta_geographic_coverage.csv")

cat("\n=== Total CTA Queries:", nrow(cov_df), "===\n")
print(
  cov_df %>%
    count(coverage_status) %>%
    mutate(pct = round(n / sum(n) * 100, 2))
)

cat("\n=== Zero-Crime Square Queries (crimes_2024 == 0) by Coverage Status ===\n")
print(
  cov_df %>%
    filter(zero_crime_query) %>%
    count(coverage_status) %>%
    mutate(pct = round(n / sum(n) * 100, 2))
)

# 2. Save visualization deliverable
plot_df <- cov_df %>%
  filter(!is.na(center_lon), !is.na(center_lat))

p <- ggplot() +
  geom_sf(data = city_metric, fill = "grey96", color = "grey45", linewidth = 0.5) +
  geom_point(
    data = plot_df,
    aes(
      x = center_lon,
      y = center_lat,
      color = coverage_status,
      shape = zero_crime_query
    ),
    size = 2.1,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "fully covered" = "#1b9e77",
      "partially covered" = "#d95f02",
      "outside Chicago" = "#7570b3"
    ),
    name = "Geographic Coverage"
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 4),
    name = "Zero-Crime Square Query\n(crimes_2024 == 0)"
  ) +
  coord_sf(crs = 4326) +
  labs(
    title = "CTA-Centered Square Query Coverage Relative to Chicago Boundary",
    subtitle = "426 CTA square queries evaluated in EPSG:26916; zero-crime flag uses crimes_2024 == 0",
    x = "Longitude",
    y = "Latitude",
    caption = "Source: City of Chicago city boundary GeoJSON & Step 4B spatial query workloads (2024)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave("figs/geographic_coverage_validation.png", plot = p, width = 9, height = 7, dpi = 300)

# 3. Save summary markdown deliverable
cov_n <- cov_df %>% count(coverage_status) %>% deframe()
z_n <- cov_df %>% filter(zero_crime_query) %>% count(coverage_status) %>% deframe()

get_cnt <- function(vec, nm) {
  val <- vec[nm]
  if (length(val) == 0 || is.na(val)) 0L else as.integer(val)
}

n_total <- nrow(cov_df)
n_full <- get_cnt(cov_n, "fully covered")
n_part <- get_cnt(cov_n, "partially covered")
n_out <- get_cnt(cov_n, "outside Chicago")

n_zero_total <- sum(cov_df$zero_crime_query, na.rm = TRUE)
n_zero_full <- get_cnt(z_n, "fully covered")
n_zero_part <- get_cnt(z_n, "partially covered")
n_zero_out <- get_cnt(z_n, "outside Chicago")

summary_md <- c(
  "# Task A: Geographic Coverage Validation Summary",
  "",
  paste0("- **Total CTA square queries evaluated:** ", n_total),
  paste0("- **Fully covered queries:** ", n_full, " (", round(n_full / n_total * 100, 2), "%)"),
  paste0("- **Partially covered queries:** ", n_part, " (", round(n_part / n_total * 100, 2), "%)"),
  paste0("- **Outside Chicago queries:** ", n_out, " (", round(n_out / n_total * 100, 2), "%)"),
  "",
  "## Zero-Crime Square Queries (`crimes_2024 == 0`)",
  "",
  paste0("- **Total zero-crime square queries:** ", n_zero_total),
  paste0("- **Outside Chicago:** ", n_zero_out),
  paste0("- **Partially covered:** ", n_zero_part),
  paste0("- **Fully covered:** ", n_zero_full),
  "",
  "## Key Takeaway",
  "",
  paste0(
    "All ", n_zero_total, " zero-crime square queries have incomplete geographic coverage ",
    "(", n_zero_out, " outside Chicago and ", n_zero_part, " partially covered), ",
    "indicating that jurisdictional boundaries likely contribute substantially to the ",
    "apparent prediction errors. For queries completely outside Chicago, the absence of ",
    "Chicago police records is expected; for partially covered queries, incomplete overlap ",
    "indicates boundary exposure rather than proving clipping alone caused a zero count."
  )
)

writeLines(summary_md, "outputs/task_a_summary.md")