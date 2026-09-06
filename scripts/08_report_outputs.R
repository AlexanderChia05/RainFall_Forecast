required_objects <- c("rain", "train", "fit", "fit_tbats", "fc", "fc_tbats", "h")
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

save_plot_quietly <- function(filename, plot, width = 8, height = 5, dpi = 150) {
  suppressMessages(
    suppressWarnings(
      ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
    )
  )
}

plot_from  <- yearmonth("2023 Jan")
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
  labs(title = "ARIMA+Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL)
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
  lo80  = as.numeric(fc_tbats$lower[, "80%"]),
  hi80  = as.numeric(fc_tbats$upper[, "80%"]),
  lo95  = as.numeric(fc_tbats$lower[, "95%"]),
  hi95  = as.numeric(fc_tbats$upper[, "95%"])
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

# Model-specific displays ---------------------------------------------------
# ETS: estimated level and seasonal states.
p_ets_states <- fit |>
  select(ets_additive) |>
  components() |>
  autoplot() +
  labs(title = "ETS(A,N,A): Estimated states")
save_plot_quietly(
  "output/plots/group_summary/ets_states.png",
  p_ets_states,
  height = 7
)

# ARIMA + Fourier: differenced series followed by its ACF and PACF.
arima_d_display <- train |>
  features(precip, unitroot_ndiffs) |>
  pull(ndiffs)

arima_series_display <- if (arima_d_display > 0) {
  train |> mutate(precip_diff = difference(precip, differences = arima_d_display))
} else {
  train |> mutate(precip_diff = precip)
}

p_diff_series <- arima_series_display |>
  autoplot(precip_diff) +
  labs(
    title = paste0("Differenced precipitation series (d = ", arima_d_display, ")"),
    x = NULL,
    y = "Differenced mm/day"
  )

p_acf_diff <- arima_series_display |>
  ACF(precip_diff, lag_max = 36) |>
  autoplot() +
  labs(title = paste0("ACF after differencing (d = ", arima_d_display, ")"))

p_pacf_diff <- arima_series_display |>
  PACF(precip_diff, lag_max = 36) |>
  autoplot() +
  labs(title = paste0("PACF after differencing (d = ", arima_d_display, ")"))

p_arima_diagnostics <- patchwork::wrap_plots(
  p_diff_series,
  p_acf_diff,
  p_pacf_diff,
  design = "AA\nBC"
) +
  patchwork::plot_annotation(
    title = "ARIMA + Fourier: Differencing and correlation diagnostics"
  )
save_plot_quietly(
  "output/plots/group_summary/arima_differenced_acf_pacf.png",
  p_arima_diagnostics,
  width = 10,
  height = 7
)

# TSLM + Fourier: regression fit and residuals versus fitted values.
tslm_aug <- augment(fit) |>
  filter(.model == "tslm_fourier") |>
  as_tibble()

p_tslm_fit <- ggplot(tslm_aug, aes(x = month)) +
  geom_line(aes(y = precip, colour = "Observed"), linewidth = 0.45) +
  geom_line(aes(y = .fitted, colour = "TSLM fitted"), linewidth = 0.55) +
  scale_colour_manual(values = c("Observed" = "grey35", "TSLM fitted" = "#0072B2")) +
  labs(title = "Observed and fitted precipitation", x = NULL, y = "mm/day", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_tslm_resid_fitted <- ggplot(tslm_aug, aes(x = .fitted, y = .innov)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.55, colour = "#D55E00", size = 1.4) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#0072B2", linewidth = 0.7) +
  labs(title = "Residuals versus fitted values", x = "Fitted mm/day", y = "Innovation") +
  theme_minimal()

p_tslm_diagnostics <- patchwork::wrap_plots(
  p_tslm_fit,
  p_tslm_resid_fitted,
  ncol = 1,
  heights = c(1.25, 1)
) +
  patchwork::plot_annotation(title = "TSLM + Fourier: Regression diagnostics")
save_plot_quietly(
  "output/plots/group_summary/tslm_fit_diagnostics.png",
  p_tslm_diagnostics,
  height = 7
)

# TBATS: estimated level and trigonometric seasonal components.
tbats_components <- forecast::tbats.components(fit_tbats)
tbats_component_names <- colnames(tbats_components)
if (is.null(tbats_component_names)) {
  tbats_component_names <- paste0("Component ", seq_len(ncol(tbats_components)))
}
tbats_components_tbl <- tibble(
  month = rep(train$month, times = ncol(tbats_components)),
  component = rep(tbats_component_names, each = nrow(tbats_components)),
  value = as.numeric(tbats_components)
)

p_tbats_components <- ggplot(tbats_components_tbl, aes(x = month, y = value)) +
  geom_line(colour = "#0072B2", linewidth = 0.45) +
  facet_wrap(vars(component), ncol = 1, scales = "free_y") +
  labs(
    title = "TBATS: Estimated model components",
    subtitle = "Components are shown on the fitted Box-Cox scale",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()
save_plot_quietly(
  "output/plots/group_summary/tbats_components.png",
  p_tbats_components,
  width = 9,
  height = 7
)

# Comparable residual panels for all four models.
p_resid_ets <- fit |>
  select(ets_additive) |>
  gg_tsresiduals()
p_resid_ets[[1]] <- p_resid_ets[[1]] +
  labs(title = "ETS(A,N,A): Residual diagnostics")
save_plot_quietly("output/plots/group_summary/resid_ets.png", p_resid_ets, height = 6)

p_resid_arima <- fit |>
  select(arima_fourier) |>
  gg_tsresiduals()
p_resid_arima[[1]] <- p_resid_arima[[1]] +
  labs(title = "ARIMA + Fourier: Residual diagnostics")
save_plot_quietly("output/plots/group_summary/resid_arima.png", p_resid_arima, height = 6)

p_resid_tslm <- fit |>
  select(tslm_fourier) |>
  gg_tsresiduals()
p_resid_tslm[[1]] <- p_resid_tslm[[1]] +
  labs(title = "TSLM + Fourier: Residual diagnostics")
save_plot_quietly("output/plots/group_summary/resid_tslm.png", p_resid_tslm, height = 6)

png("output/plots/group_summary/resid_tbats.png", width = 800, height = 600, res = 150)
tryCatch(
  forecast::checkresiduals(fit_tbats),
  finally = dev.off()
)

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
save_plot_quietly("output/plots/group_summary/fc_all_combined.png", p_all)

cat("\nDone. Wrote forecast, model-specific and residual diagnostic plots.\n")
