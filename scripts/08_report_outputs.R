dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)

plot_from  <- yearmonth("2013 Jan")
test_start <- max(train$month) + 1
test_actual_tbl <- rain |> as_tibble() |> filter(month >= test_start) |>
  transmute(month, precip)

overlay_red_actual <- function(member_name) {
  fc |> filter(.model == member_name) |>
    autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
    geom_line(data = test_actual_tbl, aes(x = month, y = precip),
              color = "red", linewidth = 0.45) +
    theme_minimal()
}

p_A <- overlay_red_actual("ets_A") +
  labs(title = "ETS: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_A_ets.png", p_A, width = 8, height = 5, dpi = 150)

p_B <- overlay_red_actual("arima_B") +
  labs(title = "ARIMA+Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_B_arima.png", p_B, width = 8, height = 5, dpi = 150)

p_C <- overlay_red_actual("tslm_C") +
  labs(title = "TSLM: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_C_tslm.png", p_C, width = 8, height = 5, dpi = 150)

tbats_train_hist <- rain |> as_tibble() |>
  filter(month >= plot_from, month < test_start) |>
  transmute(month, value = precip)
tbats_fc <- tibble(
  month = rain$month[(nrow(rain) - h + 1):nrow(rain)],
  value = as.numeric(fc_tbats$mean),
  lo80  = fc_tbats$lower[, "80%"], hi80 = fc_tbats$upper[, "80%"],
  lo95  = fc_tbats$lower[, "95%"], hi95 = fc_tbats$upper[, "95%"]
)
p_D <- ggplot() +
  geom_ribbon(data = tbats_fc, aes(month, ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(data = tbats_fc, aes(month, ymin = lo80, ymax = hi80), fill = "steelblue", alpha = 0.35) +
  geom_line(data = tbats_train_hist, aes(month, value), color = "black", linewidth = 0.6) +
  geom_line(data = tbats_fc, aes(month, value), color = "steelblue4", linewidth = 0.7) +
  geom_line(data = test_actual_tbl, aes(month, precip), color = "red", linewidth = 0.45) +
  labs(title = "TBATS: Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_D_tbats.png", p_D, width = 8, height = 5, dpi = 150)

combined_fc <- bind_rows(
  fc |> filter(.model != "snaive") |> as_tibble() |> transmute(month, .model, value = .mean),
  tibble(month = tbats_fc$month, .model = "tbats_D", value = tbats_fc$value)
)
p_all <- ggplot() +
  geom_line(data = rain |> as_tibble() |> filter(month >= plot_from), aes(month, precip),
            color = "black", linewidth = 0.6) +
  geom_line(data = combined_fc, aes(month, value, color = .model), linewidth = 0.8) +
  labs(title = "All 4 Final Models: Forecast vs Actual (test window)",
       y = "mm/day", x = NULL, color = "Model") +
  theme_minimal()
ggsave("output/plots/group_summary/fc_all_combined.png", p_all, width = 10, height = 6, dpi = 150)
