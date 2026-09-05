# 04_StephQF_arima.R - STANDALONE script, Member B (Stephanie Lo Qian
# Fui). Topic: KL monthly mean precipitation rate (mm/day), NASA POWER,
# 1981-2025 (540 obs). SDG 13 primary. Model family: dynamic-regression
# ARIMA with deterministic Fourier seasonal terms, variant
# ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0,0,0)).
#
# Self-contained: pulls its own data, runs its own diagnostics, fits and
# validates its own model, writes its own outputs. No source(), no
# readRDS of a shared file, no dependency on any other script. The setup
# and data blocks are duplicated across the four member scripts on
# purpose so each one can be submitted and run on its own.
#
# Model pick: the seasonal ARIMA search is switched OFF with PDQ(0,0,0)
# and the annual cycle is carried by 4 Fourier pairs instead. Leaving
# both in would make the seasonal differencing and the Fourier terms
# compete for the same structure. K = 4 was chosen over K = 6 because
# the extra harmonics improved in-sample fit but generalised worse.

# ---------------------------------------------------------------- setup
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

# ----------------------------------------------------------------- data
# Monthly corrected total precipitation (PRECTOTCORR, mm/day) for Kuala
# Lumpur (3.1390 N, 101.6869 E) from the NASA POWER API.
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

# ------------------------------------------------------------------ EDA
rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
rain |> gg_season(precip)    + labs(title = "Seasonal plot - monsoon cycle")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip,  lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# Trend vs seasonal strength.
print(rain |> features(precip, feat_stl))

# ------------------------------------------- stationarity / white noise
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

# --------------------------------------------------- train / test split
# Chronological hold-out, no random split: last 12 months are the test.
h     <- 12
train <- rain |> filter(month <= max(month) - h)

# --------------------------------------------------------------- model
fit <- train |> model(
  snaive      = SNAIVE(precip),
  arima_four4 = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0))
)
fc <- fit |> forecast(h = h)

print(fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE))
fit |> select(arima_four4) |> report()
fit |> select(arima_four4) |> gg_tsresiduals()

# ------------------------------------------------- residual diagnostics
# Ljung-Box degrees of freedom must be reduced by the number of AR/MA
# parameters the model itself estimated, or the test overstates how
# random the residuals are. Counted from the fitted coefficient names
# (ar1, ar2, ..., ma1, ma2, ...) rather than assumed, since PDQ(0,0,0)
# means there is no seasonal ar/ma term to also count, and the Fourier
# regressor coefficients are named after the regressor, not ar/ma, so
# this pattern only ever matches genuine ARMA error terms.
arima_coefs  <- fit |> select(arima_four4) |> tidy()
arima_dof    <- sum(grepl("^(ar|ma)[0-9]+$", arima_coefs$term))
cat("\nARIMA AR/MA parameter count used as Ljung-Box dof:", arima_dof, "\n")

cat("\n== Ljung-Box on ARIMA residuals (want p > 0.05) ==\n")
print(augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 12, dof = arima_dof))
print(augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 24, dof = arima_dof))

cat("\n== Residual ACF lags outside the white-noise band ==\n")
print(augment(fit) |> filter(.model == "arima_four4") |> as_tibble() |>
        summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
                  n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)))

# ------------------------------------- overfitting check: single holdout
acc_train <- fit |> accuracy() |> filter(.model == "arima_four4") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test  <- fc |> accuracy(rain) |> filter(.model == "arima_four4") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
holdout <- acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio_holdout = RMSE_test / RMSE_train,
         gap_pct_holdout    = abs(MASE_test - MASE_train) / MASE_test)
print(holdout)

# --------------------------------- overfitting check: rolling-origin CV
# Folds run only over the training window (up to 2024-12), never the
# 2025 holdout - otherwise the same 2025 observations that back the
# holdout accuracy above would also leak into the CV folds and inflate
# both numbers on the same data. Fit on the first 360 months, then
# refit every 6 months, each fit scored on the following 12. The CV
# ratio and gap are the authoritative overfitting numbers - the single
# holdout above is one draw only and can flatter or punish a model by
# luck of the window.
n_folds_clean <- length(seq(360, nrow(train) - h, by = 6))

set.seed(2026)
cv_acc <- train |>
  stretch_tsibble(.init = 360, .step = 6) |>
  model(arima_four4 = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0))) |>
  forecast(h = h) |>
  accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  filter(!is.na(MASE), .id <= n_folds_clean) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  member        = "arima_B",
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

# ------------------------------------------------------------- outputs
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/member_B_arima_results.csv", row.names = FALSE)

# Forecast vs actual: black history, blue forecast + interval, red actual
# overlaid on the test window so band coverage is readable.
plot_from       <- yearmonth("2013 Jan")
test_actual_tbl <- rain |> as_tibble() |> filter(month > max(train$month)) |>
  transmute(month, precip)

p_fc <- fc |> filter(.model == "arima_four4") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  geom_line(data = test_actual_tbl, aes(x = month, y = precip),
            color = "red", linewidth = 0.45) +
  labs(title = "ARIMA + Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_B_arima.png", p_fc, width = 8, height = 5, dpi = 150)

ggsave("output/plots/group_summary/resid_B_arima.png",
       fit |> select(arima_four4) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)

cat("\nDone. Wrote output/member_B_arima_results.csv and 2 plots.\n")
