# 05_HamGQ_tslm.R - STANDALONE script, Member C (Ham Guan Quan).
# Topic: KL monthly mean precipitation rate (mm/day), NASA POWER,
# 1981-2025 (540 obs). SDG 13 primary. Model family: TSLM (time series
# linear regression), variant TSLM(precip ~ trend() + fourier(K = 4)).
#
# Self-contained: pulls its own data, runs its own diagnostics, fits and
# validates its own model, writes its own outputs. No source(), no
# readRDS of a shared file, no dependency on any other script. The setup
# and data blocks are duplicated across the four member scripts on
# purpose so each one can be submitted and run on its own.
#
# Model pick: a deterministic regression on a linear trend plus 4 Fourier
# pairs. Unlike the other members this model KEEPS an explicit trend
# term, which gives a direct significance test of the long-run trend -
# its coefficient can be read against the Mann-Kendall result below. The
# two should agree: a statistically significant but practically small
# upward movement.

# setup
pkgs <- c("fpp3", "tseries", "zoo", "Kendall")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
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

# Trend vs seasonal strength - compare against the trend() coefficient.
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
    "small; 0.1-0.3 = weak. If this agrees with the trend() coefficient\n",
    "reported below (significant p, small slope), the honest reading is a\n",
    "real but practically minor upward trend, not an absence of trend.\n")

cat("\n== Min precip value (near-zero check for MAPE stability) ==\n")
print(min(rain$precip, na.rm = TRUE))

# train / test split
h     <- 12
train <- rain |> filter(month <= max(month) - h)

# model
fit <- train |> model(
  snaive = SNAIVE(precip),
  tslm   = TSLM(precip ~ trend() + fourier(K = 4))
)
fc <- fit |> forecast(h = h)

print(fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE))

fit |> select(tslm) |> report()
fit |> select(tslm) |> gg_tsresiduals()

# residual diagnostics
cat("\n== Ljung-Box on TSLM residuals (want p > 0.05) ==\n")
print(augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 12))
print(augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 24))

cat("\n== Residual ACF lags outside the white-noise band ==\n")
print(augment(fit) |> filter(.model == "tslm") |> as_tibble() |>
        summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
                  n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)))

# overfitting check: single holdout
acc_train <- fit |> accuracy() |> filter(.model == "tslm") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test  <- fc |> accuracy(rain) |> filter(.model == "tslm") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
holdout <- acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio_holdout = RMSE_test / RMSE_train,
         gap_pct_holdout    = abs(MASE_test - MASE_train) / MASE_test)
print(holdout)

# overfitting check: rolling-origin CV
n_folds_clean <- length(seq(360, nrow(rain) - h, by = 6))

set.seed(2026)
cv_acc <- rain |>
  stretch_tsibble(.init = 360, .step = 6) |>
  model(tslm = TSLM(precip ~ trend() + fourier(K = 4))) |>
  forecast(h = h) |>
  accuracy(rain, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  filter(!is.na(MASE), .id <= n_folds_clean) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  member        = "tslm_C",
  MASE_train    = holdout$MASE_train,
  RMSE_train    = holdout$RMSE_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MASE_cv    = cv_summary$sd_MASE,
  n_folds       = cv_summary$n_folds,
  rmse_ratio_cv = cv_summary$mean_RMSE / holdout$RMSE_train,
  gap_pct_cv    = abs(cv_summary$mean_MASE - holdout$MASE_train) / cv_summary$mean_MASE
) |>
  mutate(within_1_3x_cv  = rmse_ratio_cv <= 1.3,
         within_10pct_cv = gap_pct_cv <= 0.10)
print(results)

# outputs
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/member_C_tslm_results.csv", row.names = FALSE)

# Forecast vs actual
plot_from       <- yearmonth("2013 Jan")
test_actual_tbl <- rain |> as_tibble() |> filter(month > max(train$month)) |>
  transmute(month, precip)

p_fc <- fc |> filter(.model == "tslm") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  geom_line(data = test_actual_tbl, aes(x = month, y = precip),
            color = "red", linewidth = 0.45) +
  labs(title = "TSLM (trend + Fourier K=4): Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_C_tslm.png", p_fc, width = 8, height = 5, dpi = 150)

ggsave("output/plots/group_summary/resid_C_tslm.png",
       fit |> select(tslm) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)

cat("\nDone. Wrote output/member_C_tslm_results.csv and 2 plots.\n")
