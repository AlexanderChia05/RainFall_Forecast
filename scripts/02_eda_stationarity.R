source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

# EDA
rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
rain |> gg_season(precip) + labs(title = "Seasonal plot - by year")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip, lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# Stationarity + white-noise check
cat("\n== ADF (want p < 0.05 for stationary) ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS (want p > 0.05 for stationary - i.e. NOT rejecting H0) ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box on RAW series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

# Trend detection - Mann-Kendall test 
cat("\n== Mann-Kendall trend test (H0: no monotonic trend) ==\n")
print(Kendall::MannKendall(rain$precip))
cat("tau near 0 (|tau|<0.1) = negligible trend magnitude even if p is small;\n",
    "0.1-0.3 = weak; 0.3-0.5 = moderate; >0.5 = strong. Compare this tau\n",
    "against 05_HamGQ_tslm.R's trend() coefficient/p-value - if BOTH agree\n",
    "(significant p, small tau/slope), the honest conclusion is 'a real but\n",
    "practically minor upward trend', not 'no trend'.\n")

# Zero value check
cat("\n== Min precip value (near-zero check) ==\n")
print(min(rain$precip, na.rm = TRUE))
