# 04_StephQF_arima.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 primary. Member B
# model family: ARIMA with deterministic Fourier seasonal terms.
#
# Model pick: ARIMA(precip ~ fourier(K=4) + pdq() + PDQ(0,0,0)).
# PDQ(0,0,0) explicitly turns off ARIMA's own seasonal search - without
# it, fourier() and ARIMA's auto (P,D,Q)[12] search compete for the same
# 12-month signal and the fit becomes exactly singular (hit this bug in
# the trade-topic grid search, fixed proactively here from the start).
#
# K=4, not K=2 or K=6 or plain auto: tested all four (grid search since
# removed once the search phase was done - both variants' full numbers
# held here for the record):
#   arima_fourierK4 (picked): MASE_train=0.738, RMSE_train=2.44,
#     mean_MASE_cv=0.797, mean_RMSE_cv=2.51, ratio=1.03, gap=7.3%,
#     p12=0.553, p24=0.137 (final numbers per the 29-fold-clean CV -
#     an earlier pre-fold-parity-fix run showed 0.801/2.52/7.8%, same
#     ratio; see 07_group_comparison.R for the fold-count fix)
#   arima_fourierK6: MASE_train=0.737, RMSE_train=2.43,
#     mean_MASE_cv=0.804, mean_RMSE_cv=2.53, ratio=1.04, gap=8.3%,
#     p12=0.573, p24=0.122
# K=6 has the marginally better single-holdout MASE (0.814 vs K=4's
# 0.843) AND marginally better LB p-values, but K=4 has the lower
# CV-based RMSE ratio (1.03 vs 1.04) and lower MASE gap (7.8% vs 8.3%)
# - fewer parameters generalise very slightly better on the
# authoritative CV metric, same "simpler wins on CV even when a
# single-shot criterion looks better elsewhere" pattern as the ozone
# project's member B (K=1 over K=2 there). Plain arima_auto (no fourier
# at all) failed Ljung-Box outright (p=0.266 @lag12 but p=0.000632
# @lag24) - the deterministic seasonal term is doing real work here,
# not redundant with what auto-ARIMA finds alone.
#
# Diagnostics (single holdout): p=0.553 (lag12), p=0.137 (lag24) - both
# clear 0.05 with real margin. 29-fold rolling CV: RMSE ratio = 1.03 (an
# almost dead heat with 1.0, best of the 3 fable members - TBATS still
# ahead overall, see 06_ChiaZY_tbats.R), MASE gap = 7.3% - one of only 3
# picks across 3 topics this project attempted (ozone, trade, rain)
# that clears the group's original 10% overfitting target outright.

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
# 07_group_comparison.R for the authoritative CV-based ratio/gap)
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
