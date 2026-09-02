# 06_ChiaZY_tbats.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 primary. Member D
# model family: TBATS (Box-Cox, ARMA errors, Trend, Seasonal -
# trigonometric/harmonic seasonal representation, De Livera, Hyndman &
# Snyder, 2011).
#
# Replaces STL(robust)+ARIMA(remainder) (the original member D pick,
# see git history) - STL passed both LB lags (p=0.312/0.139) but had
# the weakest CV numbers of the whole final lineup (ratio=1.18,
# gap=22.7%, the only one of the 4 that missed the group's own 10%
# target). TBATS was tested afterward as a supplementary family
# (10_tbats_nnar_grid.R, since removed - grid search is done, evidence
# kept here) and outperformed every other family tested, STL included.
#
# Model pick: tbats(seasonal.periods=12, use.trend=FALSE). Verified via
# the (now-removed) TBATS/NNAR grid search - 4 variants tested:
#   tbats_notrend      : ratio=1.01, gap=5.7%, p12=0.524, p24=0.0886  <- picked
#   tbats_auto         : ratio=1.01, gap=6.4%, p12=0.539, p24=0.0994
#   tbats_boxcox       : ratio=1.02, gap=6.7%, p12=0.539, p24=0.0994
#   tbats_trend_damped : ratio=1.03, gap=7.5%, p12=0.475, p24=0.0835
# tbats_notrend has the best ratio AND lowest gap% of all 4 - matches
# ETS's own finding (gamma≈0.0001, i.e. no meaningful trend state
# needed - see 02_eda_stationarity.R's weak trend_strength=0.181):
# forcing trend off outperforms letting TBATS's own AIC-based auto-
# selection decide, same "simpler wins on CV even when a single
# in-sample criterion prefers otherwise" pattern seen everywhere else
# in this project (member B's K=4 over K=6, TSLM's fourier pick, etc).
#
# Trade-off disclosed, not hidden: TBATS's Ljung-Box lag=24 margin
# (p=0.0886) is thinner than the other 3 members (all >0.12) - it
# clears 0.05 but with less room. This is the one place TBATS is
# WORSE than what it replaced (STL had p24=0.139) - traded a thinner
# p-value margin for a much better overfitting profile (ratio 1.01 vs
# STL's 1.18, gap 5.7% vs STL's 22.7%). Flag both sides if asked to
# justify robustness in the report.
#
# Not a fable/tidyverts model - uses forecast::tbats() on a plain ts
# object, bridged in/out of the tsibble pipeline manually (see
# 00_setup.R for why forecast:: is called explicitly instead of
# library(forecast)).

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

train_ts <- ts(train$precip, frequency = 12)

fit_tbats <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                              use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats  <- forecast::forecast(fit_tbats, h = h)
print(fit_tbats)

test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats <- forecast::accuracy(fc_tbats, test_actual)
print(acc_tbats)

resid_tbats <- residuals(fit_tbats)

# Ljung-Box - residuals must look random (p > 0.05), both conventions
Box.test(resid_tbats, lag = 12, type = "Ljung-Box")
Box.test(resid_tbats, lag = 24, type = "Ljung-Box")

# ACF-in-bounds check (MUST)
acf_out_of_bounds(resid_tbats, lag.max = 12)
acf_out_of_bounds(resid_tbats, lag.max = 24)

# Overfitting check (single holdout, reference only - see
# 07_group_comparison.R for the authoritative CV-based ratio/gap)
mase_train <- acc_tbats["Training set", "MASE"]
mase_test  <- acc_tbats["Test set", "MASE"]
rmse_train <- acc_tbats["Training set", "RMSE"]
rmse_test  <- acc_tbats["Test set", "RMSE"]
cat("rmse_ratio_holdout:", rmse_test / rmse_train,
    " mase_gap_pct_holdout:", abs(mase_test - mase_train) / mase_test, "\n")

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
png("output/plots/group_summary/resid_D_tbats.png", width = 800, height = 600, res = 150)
forecast::checkresiduals(fit_tbats)
dev.off()
