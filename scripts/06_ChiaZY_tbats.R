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
# Model pick: tbats(seasonal.periods = 12, use.box.cox = FALSE,
# use.trend = FALSE). Confirmed by a 42-configuration grid search over
# Box-Cox / trend / damped-trend / arma-errors (see 09_tbats_grid.R,
# which collapses to 7 distinct fitted models). Screened on a Ljung-Box
# test whose degrees of freedom are corrected for the parameters TBATS
# itself estimated (forecast::modeldf()) - Box.test()'s fitdf defaults
# to 0, which understates how many parameters were fitted and overstates
# how random the residuals look. Under that corrected test, diagnostic
# adequacy fell as model complexity rose: every configuration with a
# Box-Cox transform, a damped trend, or both failed at lag 24 (p as low
# as 0.0124), and the single survivor was the simplest configuration in
# the family - no Box-Cox, no trend, 3 harmonics - at p = 0.0587, a thin
# but real pass. The originally-fitted configuration
# (use.box.cox = NULL, letting TBATS estimate lambda) converged to
# lambda = 0.596 and FAILED the corrected test at lag 24 (p = 0.0276),
# despite having the better test-set MASE (0.685 vs 0.730 for the
# no-Box-Cox survivor). The survivor was selected anyway: a model whose
# residuals are not white cannot support valid prediction intervals,
# so the accuracy loss was accepted in exchange for a passing
# diagnostic. Every configuration in the grid, this one included,
# converged to {0,0} AR/MA orders - TBATS never selected an ARMA error
# component on this series, with use.arma.errors TRUE or FALSE making no
# difference to the fitted model.
#
# Trade-off disclosed, not hidden: even the survivor's lag-24 margin
# (p=0.0587) is thin, and without the dof correction every configuration
# in the grid would have passed (min p=0.0835) - the pass/fail verdict
# depends entirely on which Ljung-Box convention is used, and that fork
# is reported rather than hidden behind a single number.
#
# Not a fable/tidyverts model - uses forecast::tbats() on a plain ts
# object, bridged in and out of the tsibble pipeline by hand.

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
  as_tsibble(index = month)

missing_before <- sum(is.na(rain$precip))

rain <- rain |>
  mutate(precip = zoo::na.approx(precip, na.rm = FALSE)) |>
  select(month, precip)

missing_after <- sum(is.na(rain$precip))

cat("rain:", nrow(rain), "obs,", format(min(rain$month)), "to", format(max(rain$month)), "\n")
cat("Missing before interpolation:", missing_before,
    "| Missing after interpolation:", missing_after, "\n")

# EDA
rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
rain |> gg_season(precip)    + labs(title = "Seasonal plot - monsoon cycle")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip,  lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# Seasonal-trend decomposition: the evidence behind use.trend = FALSE.
# Seasonal strength ~0.469 against trend strength ~0.181. The seasonal
# component also widens across the record, which would ordinarily
# motivate a Box-Cox transform - but the grid search above found that
# adding one costs more Ljung-Box degrees of freedom than it buys back
# in residual whiteness, so this model does not use one (see header).
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
fit_tbats <- forecast::tbats(train_ts, use.box.cox = FALSE, use.trend = FALSE,
                             use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats  <- forecast::forecast(fit_tbats, h = h)
print(fit_tbats)
cat("\nBox-Cox lambda:", if (is.null(fit_tbats$lambda)) "none (use.box.cox = FALSE)" else fit_tbats$lambda, "\n")

test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats   <- forecast::accuracy(fc_tbats, test_actual)
print(acc_tbats)

snaive_fit <- train |> model(snaive = SNAIVE(precip))
print(snaive_fit |> forecast(h = h) |> accuracy(rain) |>
        select(.model, MASE, RMSE, MAE, MAPE))

# residual diagnostics
resid_tbats <- residuals(fit_tbats)

# Ljung-Box degrees of freedom must be reduced by the number of
# parameters TBATS itself estimated (its internal ARMA error terms plus
# any Box-Cox/seasonal state count that forecast::modeldf() reports), or
# the test overstates how random the residuals are. Box.test()'s fitdf
# defaults to 0, which is why this is set explicitly rather than left
# at the default.
tbats_dof <- forecast::modeldf(fit_tbats)
cat("\nTBATS model degrees of freedom used as Ljung-Box fitdf:", tbats_dof, "\n")

cat("\n== Ljung-Box on TBATS residuals (want p > 0.05) ==\n")
print(Box.test(resid_tbats, lag = 12, type = "Ljung-Box", fitdf = tbats_dof))
print(Box.test(resid_tbats, lag = 24, type = "Ljung-Box", fitdf = tbats_dof))

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
# Folds run only over the training window (up to 2024-12), never the
# 2025 holdout - otherwise the same 2025 observations that back the
# holdout accuracy above would also leak into the CV folds and inflate
# both numbers on the same data. TBATS is not a fable model, so the
# folds are looped by hand; every slice below is drawn from train, not
# rain, so 2025 is never touched inside this loop.
origins <- seq(360, nrow(train) - h, by = 6)

cv_tbats <- map_dfr(origins, function(i) {
  tr  <- train |> slice(1:i)
  te  <- train |> slice((i + 1):(i + h)) |> pull(precip)
  m   <- forecast::tbats(ts(tr$precip, frequency = 12), use.box.cox = FALSE,
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
  lb_pvalue_24  = Box.test(resid_tbats, lag = 24, type = "Ljung-Box", fitdf = tbats_dof)$p.value
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
