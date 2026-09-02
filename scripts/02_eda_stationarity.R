# 02_eda_stationarity.R - group leader script. EDA + stationarity/white-
# noise checks for the rain topic (v3), BEFORE committing to it or
# building a family grid search - same "verify first" ordering used
# throughout this project after the trade-topic grid search came back
# 0/19 passers.

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

# EDA
rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
rain |> gg_season(precip) + labs(title = "Seasonal plot - by year")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip, lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# STL decomposition - quick look at trend/season/remainder split
rain |> model(STL(precip, robust = TRUE)) |> components() |> autoplot()

# Stationarity + white-noise check (same convention as v1/v2 - not a
# hard requirement per CONTEXT.md's own note, but must be checked and
# reported honestly either way)
cat("\n== ADF (want p < 0.05 for stationary) ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS (want p > 0.05 for stationary - i.e. NOT rejecting H0) ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box on RAW series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

# Seasonal strength (feasts::feat_stl) - quantifies how much of the
# variance is explained by the seasonal component (0-1 scale). This
# topic's whole selling point is "clean mechanical monsoon seasonality"
# - verify that claim numerically, don't just eyeball the season plot.
cat("\n== Seasonal/trend strength (feat_stl) ==\n")
print(rain |> features(precip, feat_stl))

# Zero/near-zero value check - tropical KL shouldn't have true zero-
# rain months, but confirm rather than assume (relevant for MAPE, which
# blows up near zero).
cat("\n== Min precip value (near-zero check) ==\n")
print(min(rain$precip, na.rm = TRUE))
