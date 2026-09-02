# 06_ChiaZY_tbats.R - STANDALONE script, Member D (Chia Zhen Yang).
# Topic: KL monthly mean precipitation rate (mm/day), NASA POWER,
# 1981-2025 (540 obs). SDG 13 primary. Model family: TBATS (Box-Cox,
# ARMA errors, Trend, Seasonal - trigonometric/harmonic seasonal
# representation, De Livera, Hyndman & Snyder, 2011).
#
# Self-contained: pulls its own data, runs its own diagnostics, fits and
# validates its own model, writes its own outputs. No source(), no
# readRDS of a shared file, no dependency on any other script. The setup
# and data blocks are duplicated across the four member scripts on
# purpose so each one can be submitted and run on its own.
#
# Model pick: tbats(seasonal.periods = 12, use.trend = FALSE). Verified
# via a 4-variant grid search (grid script since removed, evidence kept
# here):
#   tbats_notrend      : ratio=1.01, gap=5.7%, p24=0.0886   <- picked
#   tbats_auto         : ratio=1.01, gap=6.4%, p24=0.0994
#   tbats_boxcox       : ratio=1.02, gap=6.7%, p24=0.0994
#   tbats_trend_damped : ratio=1.03, gap=7.5%, p24=0.0835
# tbats_notrend has the best ratio AND the lowest gap% of the four,
# matching the weak STL trend strength of 0.181 measured below: forcing
# the trend off beats letting the AIC-based auto-selection decide.
#
# Trade-off disclosed, not hidden: the Ljung-Box lag=24 margin
# (p=0.0886) is thinner than the other members. It clears 0.05 but with
# less room, traded against a much better overfitting profile.
#
# Not a fable/tidyverts model - uses forecast::tbats() on a plain ts
# object, bridged in and out of the tsibble pipeline by hand. The Box-Cox
# transform is applied and reversed INSIDE the forecast package, so
# forecasts, fitted values and residuals all come back on the original
# mm/day scale and no manual back-transform is needed.

# setup
pkgs <- c("fpp3", "tseries", "zoo", "forecast", "Kendall")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
library(purrr)
library(tidyr)
library(fpp3)
library(tseries)

# Counts residual ACF lags outside the +/- 1.96/sqrt(n) white-noise band.
acf_out_of_bounds <- function(resid, lag.max = 12) {
  r  <- na.omit(resid)
  n  <- length(r)
  ci <- 1.96 / sqrt(n)
  a  <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# data
url <- paste0(
  "https://power.larc.nasa.gov/api/temporal/monthly/point?",
  "parameters=PRECTOTCORR&community=AG",
  "&longitude=101.6869&latitude=3.1390",
  "&start=1981&end=2025&format=CSV"
)

raw_lines  <- system(paste0("curl -s ", shQuote(url)), intern = TRUE)
start_line <- which(grepl("-END HEADER-", raw_lines)) + 1
df <- read.csv(text = paste(raw_lines[start_line:length(raw_lines)], collapse = "\n"),
               stringsAsFactors = FALSE)

month_abbr <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

rain <- df |>
  select(YEAR, all_of(month_abbr)) |>
  pivot_longer(-YEAR, names_to = "month_abbr", values_to = "precip") |>
  mutate(
    precip    = na_if(precip, -999),
    month_num = match(month_abbr, month_abbr),
    month     = yearmonth(paste(YEAR, month_num, sep = "-"))
  ) |>
  arrange(month) |>
  as_tsibble(index = month) |>
  mutate(precip = zoo::na.approx(precip, na.rm = FALSE)) |>
  select(month, precip)

cat("rain:", nrow(rain), "obs,", format(min(rain$month)), "to", format(max(rain$month)),
    "| NAs remaining:", sum(is.na(rain$precip)), "\n")

# EDA
rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
rain |> gg_season(precip)    + labs(title = "Seasonal plot - monsoon cycle")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip,  lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# Seasonal-trend decomposition: the evidence behind use.trend = FALSE.
# Seasonal strength ~0.469 against trend strength ~0.181, and the
# seasonal component widens across the record, which is what motivates
# letting TBATS estimate a Box-Cox transform.
p_stl <- rain |> model(STL(precip)) |> components() |> autoplot() +
  labs(title = "STL decomposition")
print(p_stl)
print(rain |> features(precip, feat_stl))

# stationarity / white noise
cat("\n== ADF (want p < 0.05 for stationary) ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS (want p > 0.05 for stationary) ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box on RAW series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

cat("\n== Mann-Kendall trend test (H0: no monotonic trend) ==\n")
print(Kendall::MannKendall(rain$precip))
cat("tau near 0 (|tau| < 0.1) = negligible trend magnitude even if p is\n",
    "small; 0.1-0.3 = weak. A significant p with a small tau means a real\n",
    "but practically minor upward trend, not an absence of trend.\n")

cat("\n== Min precip value (near-zero check for MAPE stability) ==\n")
print(min(rain$precip, na.rm = TRUE))

# train / test split
h        <- 12
train    <- rain |> filter(month <= max(month) - h)
train_ts <- ts(train$precip, frequency = 12)

# model
fit_tbats <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                             use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats  <- forecast::forecast(fit_tbats, h = h)
print(fit_tbats)
cat("\nFitted Box-Cox lambda:", fit_tbats$lambda, "\n")

test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats   <- forecast::accuracy(fc_tbats, test_actual)
print(acc_tbats)

snaive_fit <- train |> model(snaive = SNAIVE(precip))
print(snaive_fit |> forecast(h = h) |> accuracy(rain) |>
        select(.model, MASE, RMSE, MAE, MAPE))

# residual diagnostics
resid_tbats <- residuals(fit_tbats)

cat("\n== Ljung-Box on TBATS residuals (want p > 0.05) ==\n")
print(Box.test(resid_tbats, lag = 12, type = "Ljung-Box"))
print(Box.test(resid_tbats, lag = 24, type = "Ljung-Box"))

cat("\n== Residual ACF lags outside the white-noise band ==\n")
cat("n_lags_out_12:", acf_out_of_bounds(resid_tbats, lag.max = 12),
    " n_lags_out_24:", acf_out_of_bounds(resid_tbats, lag.max = 24), "\n")

# overfitting check: single holdout
mase_train <- acc_tbats["Training set", "MASE"]
mase_test  <- acc_tbats["Test set", "MASE"]
rmse_train <- acc_tbats["Training set", "RMSE"]
rmse_test  <- acc_tbats["Test set", "RMSE"]
cat("\nrmse_ratio_holdout:", rmse_test / rmse_train,
    " mase_gap_pct_holdout:", abs(mase_test - mase_train) / mase_test, "\n")

# overfitting check: rolling-origin CV
origins <- seq(360, nrow(rain) - h, by = 6)

cv_tbats <- map_dfr(origins, function(i) {
  tr  <- rain |> slice(1:i)
  te  <- rain |> slice((i + 1):(i + h)) |> pull(precip)
  m   <- forecast::tbats(ts(tr$precip, frequency = 12), use.box.cox = NULL,
                         use.trend = FALSE, use.damped.trend = FALSE,
                         seasonal.periods = 12)
  fcv <- forecast::forecast(m, h = h)
  acc <- forecast::accuracy(fcv, te)
  tibble(MASE = acc["Test set", "MASE"], RMSE = acc["Test set", "RMSE"])
})

cv_summary <- cv_tbats |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  member        = "tbats_D",
  MASE_train    = mase_train,
  RMSE_train    = rmse_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MASE_cv    = cv_summary$sd_MASE,
  n_folds       = cv_summary$n_folds,
  rmse_ratio_cv = cv_summary$mean_RMSE / rmse_train,
  gap_pct_cv    = abs(cv_summary$mean_MASE - mase_train) / cv_summary$mean_MASE,
  lb_pvalue_24  = Box.test(resid_tbats, lag = 24, type = "Ljung-Box")$p.value
) |>
  mutate(within_1_3x_cv  = rmse_ratio_cv <= 1.3,
         within_10pct_cv = gap_pct_cv <= 0.10)
print(results)

# outputs
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/member_D_tbats_results.csv", row.names = FALSE)

# STL decomposition plot (Figure 1 of the individual report).
ggsave("output/plots/group_summary/stl_D_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)

# Forecast vs actual, built by hand because TBATS is not a fable model:
# black history, blue forecast + 80/95% ribbons, red actual on top.
# Zoomed to the last ~3 years so the 12-month test window is readable -
# a 12-year window makes a 12-month forecast impossible to judge.
zoom_from       <- yearmonth("2023 Jan")
test_start      <- max(train$month) + 1
test_actual_tbl <- rain |> as_tibble() |> filter(month >= test_start) |>
  transmute(month, precip)
hist_tbl <- rain |> as_tibble() |>
  filter(month >= zoom_from, month < test_start) |>
  transmute(month, value = precip)
fc_tbl <- tibble(
  month = rain$month[(nrow(rain) - h + 1):nrow(rain)],
  value = as.numeric(fc_tbats$mean),
  lo80  = fc_tbats$lower[, "80%"], hi80 = fc_tbats$upper[, "80%"],
  lo95  = fc_tbats$lower[, "95%"], hi95 = fc_tbats$upper[, "95%"]
)

p_fc <- ggplot() +
  geom_ribbon(data = fc_tbl, aes(month, ymin = lo95, ymax = hi95),
              fill = "steelblue", alpha = 0.2) +
  geom_ribbon(data = fc_tbl, aes(month, ymin = lo80, ymax = hi80),
              fill = "steelblue", alpha = 0.35) +
  geom_line(data = hist_tbl, aes(month, value), color = "black", linewidth = 0.6) +
  geom_line(data = fc_tbl, aes(month, value), color = "steelblue4", linewidth = 0.7) +
  geom_line(data = test_actual_tbl, aes(month, precip), color = "red", linewidth = 0.45) +
  labs(title = "TBATS: Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_D_tbats.png", p_fc, width = 8, height = 5, dpi = 150)

png("output/plots/group_summary/resid_D_tbats.png", width = 800, height = 600, res = 150)
forecast::checkresiduals(fit_tbats)
dev.off()

cat("\nDone. Wrote output/member_D_tbats_results.csv and 3 plots.\n")
