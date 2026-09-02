# 09_report_outputs.R - group leader script. Report-ready deliverables:
# residual diagnostic plots + forecast-vs-actual plots (one per member +
# one combined overlay) for all 4 final models. Run AFTER
# 08_group_comparison.R - needs fit, fc, rain, train, h already in the
# environment.
#
# Simpler than the ozone project's equivalent (07_report_outputs.R):
# all 4 families here are fable-native, so every plot comes straight
# from gg_tsresiduals()/autoplot() - no forecast::-package bridging or
# manual data-frame construction needed for a non-fable member (that
# was BATS's situation in the ozone project, not an issue here).

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)

# ---- Residual diagnostic plots (tsresiduals), one per member ----
ggsave("output/plots/group_summary/resid_A_ets.png",
       fit |> select(ets_A) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
ggsave("output/plots/group_summary/resid_B_arima.png",
       fit |> select(arima_B) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
ggsave("output/plots/group_summary/resid_C_tslm.png",
       fit |> select(tslm_C) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
ggsave("output/plots/group_summary/resid_D_stl.png",
       fit |> select(stl_D) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)

# ---- Forecast vs Actual, one plot per member ----
# 540 obs (1981-2025) is too long a window to read the 12-month test
# forecast against - zoom to the last ~12 years for legibility.
plot_from <- yearmonth("2013 Jan")

p_A <- fc |> filter(.model == "ets_A") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  labs(title = "Member A - ETS: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_A_ets.png", p_A, width = 8, height = 5, dpi = 150)

p_B <- fc |> filter(.model == "arima_B") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  labs(title = "Member B - ARIMA+Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_B_arima.png", p_B, width = 8, height = 5, dpi = 150)

p_C <- fc |> filter(.model == "tslm_C") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  labs(title = "Member C - TSLM: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_C_tslm.png", p_C, width = 8, height = 5, dpi = 150)

p_D <- fc |> filter(.model == "stl_D") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  labs(title = "Member D - STL+ARIMA: Forecast vs Actual", y = "mm/day", x = NULL)
ggsave("output/plots/group_summary/fc_D_stl.png", p_D, width = 8, height = 5, dpi = 150)

# ---- Combined overlay: all 4 forecasts vs actual on one chart ----
p_all <- fc |> filter(.model != "snaive") |>
  autoplot(rain |> filter(month >= plot_from), level = NULL) +
  labs(title = "All 4 Final Models: Forecast vs Actual (test window)",
       y = "mm/day", x = NULL, colour = "Model") +
  theme_minimal()
ggsave("output/plots/group_summary/fc_all_combined.png", p_all, width = 10, height = 6, dpi = 150)

cat("Saved 4 residual plots + 4 forecast plots + 1 combined overlay to output/plots/group_summary/\n")
