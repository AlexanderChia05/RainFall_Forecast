# 06_HamGQ_tslm.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 primary. Member C
# model family: TSLM (time series linear regression: trend + seasonal
# dummies).
#
# Model pick: TSLM(precip ~ trend() + season()). Verified via
# 03_family_grid_search.R: tslm_fourier(K=4) scored marginally better
# on CV ratio (1.05 vs trend+season's 1.06) but trend+season() is kept
# as the member C pick for interpretability - the report can cite
# individual monthly seasonal coefficients directly (e.g. "November's
# coefficient is X mm/day above baseline"), which fourier's sin/cos
# pairs can't offer directly. The numeric difference (1.05 vs 1.06) is
# within noise, not a meaningful accuracy trade-off.
#
# Why TSLM works here (it failed on BOTH prior topics this project
# tried): TSLM's OLS residuals are i.i.d.-assumed with no ARMA-error
# term, so it can never whiten residuals that have leftover
# autocorrelation the deterministic mean structure didn't capture -
# that's exactly what killed it on ozone (QBO leakage) and on trade
# (business-cycle momentum, LB stat >900, p=0 on every variant tried).
# Rainfall's seasonality is close enough to a clean, additive, fixed
# monthly pattern with a genuinely weak trend (trend_strength=0.181)
# that a plain OLS mean-structure fit leaves little autocorrelated
# signal behind - the family's structural weakness (no ARMA slot)
# simply isn't triggered by this topic's residual structure.
#
# Diagnostics (single holdout): p=0.494 (lag12), p=0.132 (lag24) -
# clears both. 29-fold rolling CV: RMSE ratio = 1.06, MASE gap = 9.3%
# (within the group's 10% target).

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

rain |> autoplot(precip)
rain |> gg_season(precip)
rain |> ACF(precip, lag_max = 36) |> autoplot()

adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")

h <- 12
train <- rain |> filter(month <= max(month) - h)

fit <- train |> model(
  snaive = SNAIVE(precip),
  tslm   = TSLM(precip ~ trend() + season())
)
fc <- fit |> forecast(h = h)
fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
fc |> autoplot(rain, level = c(80, 95))

fit |> select(tslm) |> report()
fit |> select(tslm) |> gg_tsresiduals()

# Ljung-Box - residuals must look random (p > 0.05)
augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 12)
augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check (MUST)
augment(fit) |> filter(.model == "tslm") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check (single holdout, reference only - see
# 08_group_comparison.R for the authoritative CV-based ratio/gap)
acc_train <- fit |> accuracy() |> filter(.model == "tslm") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "tslm") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/HamGQ_tslm", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/HamGQ_tslm/resid_tslm.png",
       fit |> select(tslm) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
