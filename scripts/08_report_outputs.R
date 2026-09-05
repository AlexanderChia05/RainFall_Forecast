required_objects <- c("rain", "train", "fc", "fc_tbats", "h")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0L) {
  stop(
    "Run scripts/07_group_comparison.R in this R session before ",
    "scripts/08_report_outputs.R. Missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)

plot_from  <- yearmonth("2013 Jan")
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
ggsave("output/plots/group_summary/fc_ets.png", p_ets, width = 8, height = 5, dpi = 150)

p_arima <- overlay_red_actual("arima_fourier") +
  labs(title = "ARIMA+Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_arima.png", p_arima, width = 8, height = 5, dpi = 150)

p_tslm <- overlay_red_actual("tslm_fourier") +
  labs(title = "TSLM: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_tslm.png", p_tslm, width = 8, height = 5, dpi = 150)

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
ggsave("output/plots/group_summary/fc_tbats.png", p_tbats, width = 8, height = 5, dpi = 150)

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
p_all <- ggplot(comparison_series, aes(month, value, color = series, linetype = series)) +
  geom_line(linewidth = 0.65) +
  geom_vline(xintercept = as.Date(test_start), linetype = "dashed",
             color = "grey35") +
  scale_color_manual(values = c(
    "Training data" = "#0072B2", "Actual test data" = "black",
    "ETS" = "#D55E00", "ARIMA + Fourier" = "#009E73",
    "TSLM + Fourier" = "#CC79A7", "TBATS" = "#E69F00"
  )) +
  scale_linetype_manual(values = c(
    "Training data" = "solid", "Actual test data" = "solid",
    "ETS" = "dashed", "ARIMA + Fourier" = "dotted",
    "TSLM + Fourier" = "dotdash", "TBATS" = "longdash"
  )) +
  labs(title = "Forecast comparison on the test set", y = "mm/day", x = NULL,
       color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("output/plots/group_summary/fc_all_combined.png", p_all, width = 8, height = 5, dpi = 150)
