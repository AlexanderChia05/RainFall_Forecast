# 05_StephQF_arima.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 primary. Member B
# model family: ARIMA with deterministic Fourier seasonal terms.
#
# Model pick: ARIMA(precip ~ fourier(K=4) + pdq() + PDQ(0,0,0)).
# PDQ(0,0,0) explicitly turns off ARIMA's own seasonal search - without
# it, fourier() and ARIMA's auto (P,D,Q)[12] search compete for the same
# 12-month signal and the fit becomes exactly singular (hit this bug in
# the trade-topic grid search, fixed proactively here from the start).
#
# K=4, not K=2 or K=6 or plain auto: tested all four in
# 03_family_grid_search.R. K=6 (max, Nyquist) scored marginally higher
# MASE (0.814 vs K=4's 0.843 on single holdout) but K=4 has the lower
# CV-based RMSE ratio (1.03 vs K=6's 1.04) and lower MASE gap (7.8% vs
# 8.3%) - fewer parameters generalise very slightly better, same
# "simpler wins on CV even when single-holdout MASE looks better"
# pattern as the ozone project's member B (K=1 over K=2 there). Plain
# arima_auto (no fourier at all) failed Ljung-Box outright (p=0.266
# @lag12 but p=0.000632 @lag24) - the deterministic seasonal term is
# doing real work here, not redundant with what auto-ARIMA finds alone.
#
# Diagnostics (single holdout): p=0.553 (lag12), p=0.137 (lag24) - both
# clear 0.05 with real margin. 29-fold rolling CV: RMSE ratio = 1.03 (an
# almost dead heat with 1.0, best of all 17 grid variants tested across
# all 4 families), MASE gap = 7.8% - the only pick across 3 topics this
# project attempted (ozone, trade, rain) that clears the group's
# original 10% overfitting target outright.

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

rain |> autoplot(precip)
rain |> gg_season(precip)
rain |> ACF(precip, lag_max = 36) |> autoplot()
rain |> PACF(precip, lag_max = 36) |> autoplot()

adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")

h <- 12
train <- rain |> filter(month <= max(month) - h)

fit <- train |> model(
  snaive      = SNAIVE(precip),
  arima_four4 = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0))
)
fc <- fit |> forecast(h = h)
fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
fc |> autoplot(rain, level = c(80, 95))

fit |> select(arima_four4) |> report()
fit |> select(arima_four4) |> gg_tsresiduals()

# Ljung-Box - residuals must look random (p > 0.05)
augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 12)
augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check (MUST)
augment(fit) |> filter(.model == "arima_four4") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check (single holdout, reference only - see
# 08_group_comparison.R for the authoritative CV-based ratio/gap)
acc_train <- fit |> accuracy() |> filter(.model == "arima_four4") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "arima_four4") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/StephQF_arima", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/StephQF_arima/resid_arima_four4.png",
       fit |> select(arima_four4) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
