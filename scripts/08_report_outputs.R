required_objects <- c("rain", "train", "fc", "fc_tbats", "cv_acc", "cv_tbats", "h")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0L) {
  source("scripts/07_group_comparison.R")
}
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0L) {
  stop("Unable to initialize report objects: ", paste(missing_objects, collapse = ", "))
}

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)

save_plot_quietly <- function(filename, plot, width = 8, height = 5, dpi = 150) {
  suppressMessages(
    suppressWarnings(
      ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
    )
  )
}

plot_from  <- yearmonth("2021 Jan")
test_start <- max(train$month) + 1
test_actual_tbl <- rain |> as_tibble() |> filter(month >= test_start) |>
  transmute(month, precip)

overlay_red_actual <- function(model_name) {
  fc |> filter(.model == model_name) |>
    autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
    geom_line(data = test_actual_tbl, aes(x = month, y = precip),
              color = "red", linewidth = 0.45) +
    theme_minimal()
}

p_ets <- overlay_red_actual("ets_additive") +
  labs(title = "ETS: Forecast vs Actual", y = "mm/day", x = NULL)
save_plot_quietly("output/plots/group_summary/fc_ets.png", p_ets)

p_arima <- overlay_red_actual("arima_fourier") +
  labs(title = "ARIMA+Fourier(K=3): Forecast vs Actual", y = "mm/day", x = NULL)
save_plot_quietly("output/plots/group_summary/fc_arima.png", p_arima)

p_tslm <- overlay_red_actual("tslm_fourier") +
  labs(title = "TSLM: Forecast vs Actual", y = "mm/day", x = NULL)
save_plot_quietly("output/plots/group_summary/fc_tslm.png", p_tslm)

tbats_train_hist <- rain |> as_tibble() |>
  filter(month >= plot_from, month < test_start) |>
  transmute(month, value = precip)
tbats_fc <- tibble(
  month = rain$month[(nrow(rain) - h + 1):nrow(rain)],
  value = as.numeric(fc_tbats$mean),
  lo80  = fc_tbats$lower[, "80%"], hi80 = fc_tbats$upper[, "80%"],
  lo95  = fc_tbats$lower[, "95%"], hi95 = fc_tbats$upper[, "95%"]
)
p_tbats <- ggplot() +
  geom_ribbon(data = tbats_fc, aes(month, ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(data = tbats_fc, aes(month, ymin = lo80, ymax = hi80), fill = "steelblue", alpha = 0.35) +
  geom_line(data = tbats_train_hist, aes(month, value), color = "black", linewidth = 0.6) +
  geom_line(data = tbats_fc, aes(month, value), color = "steelblue4", linewidth = 0.7) +
  geom_line(data = test_actual_tbl, aes(month, precip), color = "red", linewidth = 0.45) +
  labs(title = "TBATS: Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
save_plot_quietly("output/plots/group_summary/fc_tbats.png", p_tbats)

combined_fc <- bind_rows(
  fc |> filter(.model != "seasonal_naive") |> as_tibble() |>
    transmute(month, series = case_when(
      as.character(.model) == "ets_additive" ~ "ETS",
      as.character(.model) == "arima_fourier" ~ "ARIMA + Fourier",
      as.character(.model) == "tslm_fourier" ~ "TSLM + Fourier",
      TRUE ~ as.character(.model)
    ),
              value = .mean),
  tibble(month = tbats_fc$month, series = "TBATS", value = tbats_fc$value)
)

comparison_series <- bind_rows(
  rain |> as_tibble() |>
    filter(month >= plot_from, month < test_start) |>
    transmute(month, series = "Training data", value = precip),
  test_actual_tbl |>
    transmute(month, series = "Actual test data", value = precip),
  combined_fc
) |>
  mutate(series = factor(
    series,
    levels = c("Training data", "Actual test data", "ETS",
               "ARIMA + Fourier", "TSLM + Fourier", "TBATS")
  ))

# Combined forecast plot
p_all <- ggplot(comparison_series, aes(month, value, color = series)) +
  geom_line(linewidth = 0.65) +
  geom_vline(xintercept = as.Date(test_start), linetype = "dashed",
             color = "grey35") +
  scale_color_manual(values = c(
    "Training data" = "#0072B2", "Actual test data" = "black",
    "ETS" = "#D55E00", "ARIMA + Fourier" = "#009E73",
    "TSLM + Fourier" = "#CC79A7", "TBATS" = "#E69F00"
  )) +
  scale_x_yearmonth(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(title = "Forecast comparison on the test set", y = "mm/day", x = NULL,
       color = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot_quietly(
  "output/plots/group_summary/fc_all_combined.png", p_all,
  width = 4.0, height = 2.7, dpi = 220
)

# Group-report-specific summaries -----------------------------------------
dir.create("output/plots/group_report", recursive = TRUE, showWarnings = FALSE)

# A readable seasonal profile for IEEE single-column width.
month_levels <- month.abb
seasonal_profile <- rain |>
  as_tibble() |>
  mutate(calendar_month = factor(
    format(as.Date(month), "%b"),
    levels = month_levels
  )) |>
  group_by(calendar_month) |>
  summarise(
    mean = mean(precip),
    median = median(precip),
    q25 = quantile(precip, 0.25),
    q75 = quantile(precip, 0.75),
    .groups = "drop"
  )

p_seasonal_profile <- ggplot(
  seasonal_profile,
  aes(x = calendar_month, y = mean, group = 1)
) +
  geom_ribbon(aes(ymin = q25, ymax = q75, group = 1),
              fill = "#BFD7EA", alpha = 0.7) +
  geom_line(colour = "#1F4E78", linewidth = 0.8) +
  geom_point(colour = "#1F4E78", size = 1.8) +
  geom_line(aes(y = median), colour = "#D55E00", linewidth = 0.65,
            linetype = "dashed") +
  labs(
    title = "Monthly rainfall profile",
    subtitle = "1981-2025; mean, median, and IQR",
    x = NULL,
    y = "Precipitation (mm/day)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_plot_quietly(
  "output/plots/group_report/monthly_seasonal_profile.png",
  p_seasonal_profile,
  width = 4.0,
  height = 2.6,
  dpi = 220
)

# Fold-level errors make stability visible instead of reporting only means.
cv_fold_results <- bind_rows(
  cv_acc |>
    as_tibble() |>
    transmute(
      fold = as.integer(.id),
      model = recode(
        .model,
        seasonal_naive = "Seasonal naive",
        ets_additive = "ETS",
        arima_fourier = "ARIMA + Fourier",
        tslm_fourier = "TSLM + Fourier"
      ),
      MASE
    ),
  cv_tbats |>
    mutate(fold = row_number(), model = "TBATS") |>
    select(fold, model, MASE)
) |>
  arrange(fold, model)

write.csv(cv_fold_results, "output/tables/cv_fold_results.csv", row.names = FALSE)

p_cv_stability <- cv_fold_results |>
  filter(model != "Seasonal naive") |>
  ggplot(aes(x = fold, y = MASE, colour = model)) +
  geom_hline(yintercept = 1, colour = "grey45", linetype = "dashed") +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 27),
    limits = c(1, 27)
  ) +
  scale_colour_manual(values = c(
    "ETS" = "#0072B2",
    "ARIMA + Fourier" = "#009E73",
    "TSLM + Fourier" = "#CC79A7",
    "TBATS" = "#D55E00"
  )) +
  labs(
    title = "Rolling-origin MASE across 27 folds",
    subtitle = "Below 1 improves on seasonal naive",
    x = "Forecast origin",
    y = "MASE",
    colour = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

save_plot_quietly(
  "output/plots/group_report/cv_mase_by_fold.png",
  p_cv_stability,
  width = 4.0,
  height = 2.7,
  dpi = 220
)

# Refit the selected specification to all observations and issue the forecast
# from the December 2025 origin. Bias adjustment gives the conditional mean.
full_start <- c(
  as.integer(format(min(rain$month), "%Y")),
  as.integer(format(min(rain$month), "%m"))
)
full_ts <- ts(rain$precip, start = full_start, frequency = 12)
fit_tbats_full <- forecast::tbats(
  full_ts,
  use.box.cox = NULL,
  use.trend = FALSE,
  use.damped.trend = FALSE,
  seasonal.periods = 12
)
fc_tbats_2026 <- forecast::forecast(
  fit_tbats_full,
  h = 12,
  level = c(80, 95),
  biasadj = TRUE
)

future_forecast <- tibble(
  month = as.character(max(rain$month) + seq_len(12)),
  mean_forecast = as.numeric(fc_tbats_2026$mean),
  lower_80 = as.numeric(fc_tbats_2026$lower[, "80%"]),
  upper_80 = as.numeric(fc_tbats_2026$upper[, "80%"]),
  lower_95 = as.numeric(fc_tbats_2026$lower[, "95%"]),
  upper_95 = as.numeric(fc_tbats_2026$upper[, "95%"])
)
write.csv(future_forecast, "output/tables/tbats_2026_forecast.csv", row.names = FALSE)

history_tbl <- rain |>
  as_tibble() |>
  filter(month >= yearmonth("2021 Jan")) |>
  transmute(month, precip)
future_plot_tbl <- future_forecast |> mutate(month = yearmonth(month))

p_future <- ggplot() +
  geom_ribbon(data = future_plot_tbl,
              aes(month, ymin = lower_95, ymax = upper_95),
              fill = "#9ECAE1", alpha = 0.45) +
  geom_ribbon(data = future_plot_tbl,
              aes(month, ymin = lower_80, ymax = upper_80),
              fill = "#3182BD", alpha = 0.35) +
  geom_line(data = history_tbl, aes(month, precip),
            colour = "grey30", linewidth = 0.55) +
  geom_line(data = future_plot_tbl, aes(month, mean_forecast),
            colour = "#D55E00", linewidth = 0.9) +
  geom_vline(xintercept = as.Date(yearmonth("2026 Jan")),
             colour = "grey40", linetype = "dashed") +
  labs(
    title = "TBATS forecast for 2026",
    subtitle = "Dec 2025 origin; mean with 80% and 95% intervals",
    x = NULL,
    y = "Precipitation (mm/day)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_plot_quietly(
  "output/plots/group_report/tbats_2026_forecast.png",
  p_future,
  width = 4.0,
  height = 2.7,
  dpi = 220
)

cat("\nDone. Wrote report plots, fold-level errors, and final forecast.\n")
