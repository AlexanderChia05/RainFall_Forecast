# 08_report_outputs.R - group leader script. Report-ready forecast-vs-
# actual plots for all 4 final models, in the classic textbook style:
# black = training history, blue = forecast + 80/95% interval ribbon,
# red = actual overlaid ON TOP of the forecast in the test window (so
# it's visible whether the actual falls inside/outside the band). Run
# AFTER 07_group_comparison.R - needs fit, fc, fc_tbats, rain, train, h
# already in the environment.

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)

# Zoom to the last ~12 years for legibility (540 obs since 1981 is too
# long a window to read a 12-month test forecast against).
plot_from  <- yearmonth("2013 Jan")
test_start <- max(train$month) + 1
test_actual_tbl <- rain |> as_tibble() |> filter(month >= test_start) |>
  transmute(month, precip)

# ---- A/B/C: fable's own autoplot (black history + blue forecast/ribbon
# by default) + a red line layered on top for the test-window actual ----
overlay_red_actual <- function(member_name) {
  fc |> filter(.model == member_name) |>
    autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
    geom_line(data = test_actual_tbl, aes(x = month, y = precip),
              color = "red", linewidth = 0.45) +
    theme_minimal()
}

p_A <- overlay_red_actual("ets_A") +
  labs(title = "Member A - ETS: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_A_ets.png", p_A, width = 8, height = 5, dpi = 150)

p_B <- overlay_red_actual("arima_B") +
  labs(title = "Member B - ARIMA+Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_B_arima.png", p_B, width = 8, height = 5, dpi = 150)

p_C <- overlay_red_actual("tslm_C") +
  labs(title = "Member C - TSLM: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_C_tslm.png", p_C, width = 8, height = 5, dpi = 150)

# ---- D: TBATS (forecast::, not fable) - manual data frame, same 3-
# colour convention built by hand (fc_tbats$upper/$lower give the
# ribbon; test_actual_tbl gives the red overlay). Replaces the old
# fc_D_stl.png. ----
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
  labs(title = "Member D - TBATS: Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_D_tbats.png", p_D, width = 8, height = 5, dpi = 150)

# Remove the stale STL-era file if it's still sitting in the folder from
# before the member D swap - harmless if it was never there.
if (file.exists("output/plots/group_summary/fc_D_stl.png")) {
  file.remove("output/plots/group_summary/fc_D_stl.png")
  cat("Removed stale fc_D_stl.png (replaced by fc_D_tbats.png).\n")
}

# ---- Combined overlay: all 4 forecasts vs actual on one chart (no
# ribbons here - 4 overlapping bands would be unreadable, point
# forecasts only, kept from the previous version) ----
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

cat("Saved 4 forecast-vs-actual plots (black/blue/red style) + 1 combined overlay",
    "to output/plots/group_summary/\n")
