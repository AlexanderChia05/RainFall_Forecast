# 07_ChiaZY_stl.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 primary. Member D
# model family: STL decomposition forecasting.
#
# Model pick: STL(precip, robust=TRUE) + ARIMA(remainder). Verified via
# 03_family_grid_search.R: stl_notrobust scored better on single-holdout
# MASE (0.797 vs robust's 0.811) but FAILED the CV-based check
# (rmse_ratio=1.31, just over the 1.3x line, flagged Overfitting) while
# the robust version passed cleanly (rmse_ratio=1.18). Textbook case of
# why the group standardised on CV over single-holdout: robust=TRUE
# down-weights the handful of extreme-rainfall months (tropical
# downpour outliers) when estimating trend/season instead of letting
# them distort the decomposition - exactly what it's meant to do, and
# here it measurably improves generalisation, not just single-holdout
# vanity metrics.
#
# Diagnostics (single holdout): p=0.312 (lag12), p=0.139 (lag24) -
# clears both, though the thinnest margin of the 4 final picks (still
# far more comfortable than anything achieved on the ozone or trade
# topics). 29-fold rolling CV: RMSE ratio = 1.18, MASE gap = 22.7% -
# the only one of the 4 final picks that does NOT clear the group's own
# 10% target, though still comfortably within the 1.3x hard rule.
# Report honestly: STL's remainder-ARIMA step evidently captures the
# monsoon cycle less tightly than the other 3 families' more direct
# approaches, even though it's the family best equipped to shrug off
# outlier months in the estimation stage.

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

rain |> autoplot(precip)
rain |> gg_season(precip)
rain |> model(STL(precip, robust = TRUE)) |> components() |> autoplot()

adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")

h <- 12
train <- rain |> filter(month <= max(month) - h)

fit <- train |> model(
  snaive    = SNAIVE(precip),
  stl_arima = decomposition_model(STL(precip, robust = TRUE),
                                   ARIMA(season_adjust ~ pdq()))
)
fc <- fit |> forecast(h = h)
fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
fc |> autoplot(rain, level = c(80, 95))

fit |> select(stl_arima) |> gg_tsresiduals()

# Ljung-Box - residuals must look random (p > 0.05)
augment(fit) |> filter(.model == "stl_arima") |> features(.innov, ljung_box, lag = 12)
augment(fit) |> filter(.model == "stl_arima") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check (MUST)
augment(fit) |> filter(.model == "stl_arima") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check (single holdout, reference only - see
# 08_group_comparison.R for the authoritative CV-based ratio/gap)
acc_train <- fit |> accuracy() |> filter(.model == "stl_arima") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "stl_arima") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/ChiaZY_stl", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/ChiaZY_stl/resid_stl_arima.png",
       fit |> select(stl_arima) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
