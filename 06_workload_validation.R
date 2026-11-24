library(tidyverse)

# 1. Load datasets ------------------------------------------------------------
w     <- read_csv("data/processed/04b_spatial_query_workloads_2024.csv", show_col_types = FALSE)
preds <- read_csv("data/processed/05_model_predictions_2024.csv", show_col_types = FALSE)
eval5 <- read_csv("outputs/05_evaluation_metrics.csv", show_col_types = FALSE)
diag5 <- read_csv("outputs/05_area_overlap_diagnostic.csv", show_col_types = FALSE)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# Standardize workload group labels using exact query_source values
w_clean <- w %>%
  mutate(
    workload_group = case_when(
      query_source == "CTA_centered_rectangle" ~ "CTA",
      query_source == "UCI_original_rectangle" ~ "UCI",
      TRUE ~ query_source
    )
  )

preds_clean <- preds %>%
  mutate(
    workload_group = case_when(
      query_source == "CTA_centered_rectangle" |
        str_detect(evaluation_set, regex("^CTA", ignore_case = TRUE)) ~ "CTA",
      query_source == "UCI_original_rectangle" |
        str_detect(evaluation_set, regex("^UCI", ignore_case = TRUE)) ~ "UCI",
      TRUE ~ query_source
    ),
    err_base = pred_baseline_count - crimes_2024,
    err_aug  = pred_augmented_count - crimes_2024,
    ae_base  = abs(err_base),
    ae_aug   = abs(err_aug),
    se_base  = err_base^2,
    se_aug   = err_aug^2,
    ae_per_km2_base = ae_base / pmax(area_km2, 1e-6),
    ae_per_km2_aug  = ae_aug  / pmax(area_km2, 1e-6)
  )

calc_metrics <- function(df, err_col, ae_col, se_col, log_ae_col) {
  err    <- df[[err_col]]
  ae     <- df[[ae_col]]
  se     <- df[[se_col]]
  log_ae <- df[[log_ae_col]]
  mean_y <- mean(df$crimes_2024, na.rm = TRUE)
  mae_val <- mean(ae, na.rm = TRUE)
  
  tibble(
    n_queries              = nrow(df),
    mean_area_km2          = mean(df$area_km2, na.rm = TRUE),
    median_area_km2        = median(df$area_km2, na.rm = TRUE),
    mean_actual_crimes     = mean_y,
    mae_count              = mae_val,
    med_ae_count           = median(ae, na.rm = TRUE),
    rmse_count             = sqrt(mean(se, na.rm = TRUE)),
    bias_mean_error        = mean(err, na.rm = TRUE),
    rel_mae_pct_mean_crime = if_else(!is.na(mean_y) & mean_y > 0, 100 * mae_val / mean_y, NA_real_),
    mae_log1p              = mean(log_ae, na.rm = TRUE)
  )
}

# 2. Section 1: Compare UCI and CTA workloads ---------------------------------
uci_w <- w_clean %>% filter(workload_group == "UCI", !is.na(area_km2))
cta_w <- w_clean %>% filter(workload_group == "CTA", !is.na(area_km2))

stopifnot(nrow(uci_w) > 0, nrow(cta_w) > 0)

ks_area  <- suppressWarnings(ks.test(uci_w$area_km2, cta_w$area_km2))
mw_area  <- suppressWarnings(wilcox.test(uci_w$area_km2, cta_w$area_km2))
mw_lat   <- suppressWarnings(wilcox.test(uci_w$center_lat, cta_w$center_lat))
mw_lon   <- suppressWarnings(wilcox.test(uci_w$center_lon, cta_w$center_lon))
mw_crime <- suppressWarnings(wilcox.test(uci_w$crimes_2024, cta_w$crimes_2024))

workload_dist_summary <- w_clean %>%
  group_by(workload_group) %>%
  summarise(
    n_queries          = n(),
    area_min           = min(area_km2, na.rm = TRUE),
    area_q25           = quantile(area_km2, 0.25, na.rm = TRUE),
    area_median        = median(area_km2, na.rm = TRUE),
    area_q75           = quantile(area_km2, 0.75, na.rm = TRUE),
    area_max           = max(area_km2, na.rm = TRUE),
    lat_median         = median(center_lat, na.rm = TRUE),
    lat_iqr            = IQR(center_lat, na.rm = TRUE),
    lon_median         = median(center_lon, na.rm = TRUE),
    lon_iqr            = IQR(center_lon, na.rm = TRUE),
    crimes_median      = median(crimes_2024, na.rm = TRUE),
    crimes_per_km2_med = median(crimes_per_km2_2024, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ks_area_D   = as.numeric(ks_area$statistic),
    ks_area_p   = ks_area$p.value,
    mw_area_p   = mw_area$p.value,
    mw_lat_p    = mw_lat$p.value,
    mw_lon_p    = mw_lon$p.value,
    mw_crimes_p = mw_crime$p.value
  )

# 3. Section 2: Model extrapolation (189 in-range vs. 237 out-of-range CTA) ---
cta_preds <- preds_clean %>% filter(workload_group == "CTA")

if (any(!is.na(cta_preds$within_uci_training_support))) {
  cta_extrap <- cta_preds %>%
    mutate(
      extrap_status = case_when(
        isTRUE(within_uci_training_support)  ~ "In_UCI_Range",
        isFALSE(within_uci_training_support) ~ "Outside_UCI_Range",
        TRUE                                 ~ "Unknown_Support"
      )
    )
} else {
  uci_ref <- w_clean %>% filter(workload_group == "UCI")
  cta_extrap <- cta_preds %>%
    mutate(
      in_range = area_km2 >= min(uci_ref$area_km2, na.rm = TRUE) &
        area_km2 <= max(uci_ref$area_km2, na.rm = TRUE) &
        center_lat >= min(uci_ref$center_lat, na.rm = TRUE) &
        center_lat <= max(uci_ref$center_lat, na.rm = TRUE) &
        center_lon >= min(uci_ref$center_lon, na.rm = TRUE) &
        center_lon <= max(uci_ref$center_lon, na.rm = TRUE),
      extrap_status = if_else(in_range, "In_UCI_Range", "Outside_UCI_Range")
    )
}

extrap_metrics <- bind_rows(
  cta_extrap %>%
    group_by(extrap_status) %>%
    group_modify(~ calc_metrics(.x, "err_base", "ae_base", "se_base", "abs_log_error_baseline")) %>%
    mutate(model = "baseline", .before = 2),
  cta_extrap %>%
    group_by(extrap_status) %>%
    group_modify(~ calc_metrics(.x, "err_aug", "ae_aug", "se_aug", "abs_log_error_augmented")) %>%
    mutate(model = "augmented", .before = 2)
) %>%
  ungroup()

# 4. Section 3: Area-matched comparison (accounting for size & sample size) ---
area_breaks <- c(-Inf, 0.5, 1.5, 3.0, 10.0, Inf)
area_labels <- c("<=0.5_km2", "0.5_to_1.5_km2", "1.5_to_3.0_km2", "3.0_to_10.0_km2", ">10.0_km2")

matched_df <- preds_clean %>%
  mutate(area_bin = cut(area_km2, breaks = area_breaks, labels = area_labels, right = TRUE))

min_sample_threshold <- 15

area_matched_metrics <- bind_rows(
  matched_df %>%
    group_by(area_bin, workload_group) %>%
    group_modify(~ calc_metrics(.x, "err_base", "ae_base", "se_base", "abs_log_error_baseline")) %>%
    mutate(model = "baseline", .before = 3),
  matched_df %>%
    group_by(area_bin, workload_group) %>%
    group_modify(~ calc_metrics(.x, "err_aug", "ae_aug", "se_aug", "abs_log_error_augmented")) %>%
    mutate(model = "augmented", .before = 3)
) %>%
  ungroup() %>%
  mutate(insufficient_sample_flag = n_queries < min_sample_threshold)

# 5. Section 4: Evaluate whether CTA features help (baseline vs. augmented) ---
cta_feature_impact <- cta_extrap %>%
  mutate(
    ae_diff  = ae_aug - ae_base,
    log_diff = abs_log_error_augmented - abs_log_error_baseline
  ) %>%
  group_by(extrap_status) %>%
  summarise(
    n_queries                = n(),
    mean_stations            = mean(stations_inside_query, na.rm = TRUE),
    mean_cta_entries         = mean(total_cta_entries_at_stations_inside_query_2024, na.rm = TRUE),
    base_mae_count           = mean(ae_base, na.rm = TRUE),
    aug_mae_count            = mean(ae_aug, na.rm = TRUE),
    mae_count_change         = aug_mae_count - base_mae_count,
    mae_count_pct_change     = 100 * mae_count_change / base_mae_count,
    base_mae_log1p           = mean(abs_log_error_baseline, na.rm = TRUE),
    aug_mae_log1p            = mean(abs_log_error_augmented, na.rm = TRUE),
    mae_log1p_change         = aug_mae_log1p - base_mae_log1p,
    pct_queries_improved_log = 100 * mean(log_diff < 0, na.rm = TRUE),
    pct_queries_improved_cnt = 100 * mean(ae_diff < 0, na.rm = TRUE),
    wilcox_paired_log_p      = suppressWarnings(
      wilcox.test(abs_log_error_augmented, abs_log_error_baseline, paired = TRUE)$p.value
    ),
    .groups = "drop"
  )

# 6. Save CSV deliverables ----------------------------------------------------
combined_metrics <- bind_rows(
  extrap_metrics %>%
    mutate(section = "2_extrapolation", subgroup = extrap_status) %>%
    select(section, subgroup, model, everything(), -extrap_status),
  area_matched_metrics %>%
    mutate(
      section = "3_area_matched",
      subgroup = paste0(area_bin, " | ", workload_group, if_else(insufficient_sample_flag, " [LOW_N]", ""))
    ) %>%
    select(section, subgroup, model, everything(), -area_bin, -workload_group)
)

write_csv(combined_metrics, "outputs/06_workload_validation_metrics.csv")
write_csv(workload_dist_summary, "outputs/06_workload_distribution_comparison.csv")
write_csv(cta_feature_impact, "outputs/06_cta_feature_impact.csv")

# 7. Visualization (PNG) ------------------------------------------------------
p_extrap <- extrap_metrics %>%
  ggplot(aes(x = extrap_status, y = mae_count, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = paste0(round(mae_count, 0), "\n(n=", n_queries, ")")),
    position = position_dodge(width = 0.75),
    vjust = -0.2,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "CTA Extrapolation & Model Impact",
    subtitle = "MAE (crime count) by UCI feature support & model",
    x = NULL,
    y = "Mean Absolute Error (Count)",
    fill = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

p_area <- area_matched_metrics %>%
  filter(model == "augmented") %>%
  ggplot(aes(x = area_bin, y = rel_mae_pct_mean_crime, fill = workload_group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = paste0("n=", n_queries, if_else(insufficient_sample_flag, "*", ""))),
    position = position_dodge(width = 0.75),
    vjust = -0.3,
    size = 2.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Size-Normalized Error Across Area Tiers",
    subtitle = "Relative MAE (% of mean crimes) | * n < 15 (low sample)",
    x = "Query Area Tier",
    y = "Relative MAE (% of Mean Crimes)",
    fill = "Workload"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

if (requireNamespace("gridExtra", quietly = TRUE)) {
  combined_plot <- gridExtra::grid.arrange(p_extrap, p_area, ncol = 2)
  ggsave("figs/06_workload_validation.png", combined_plot, width = 12, height = 5.5, dpi = 300)
} else {
  ggsave("figs/06_workload_validation.png", p_extrap, width = 10, height = 5.5, dpi = 300)
}

print(workload_dist_summary)
print(extrap_metrics)
print(area_matched_metrics, n = Inf)
print(cta_feature_impact)