# ARIMA with Fourier terms rainfall forecast

# Setup
pkgs <- c("fpp3", "tseries", "zoo", "Kendall", "patchwork")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
library(tidyr)
library(fpp3)
library(tseries)

acf_out_of_bounds <- function(resid, lag.max = 12) {
  r  <- na.omit(resid)
  n  <- length(r)
  ci <- 1.96 / sqrt(n)
  a  <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# Data
url <- paste0(
  "https://power.larc.nasa.gov/api/temporal/monthly/point?",
  "parameters=PRECTOTCORR&community=AG",
  "&longitude=101.6869&latitude=3.1390",
  "&start=1981&end=2025&format=CSV"
)

raw_lines  <- system(paste0("curl -s ", shQuote(url)), intern = TRUE)
header_end <- which(grepl("-END HEADER-", raw_lines))
if (length(header_end) != 1L) {
  stop("NASA POWER download failed or returned an unexpected response.")
}
start_line <- header_end + 1L
df <- read.csv(text = paste(raw_lines[start_line:length(raw_lines)], collapse = "\n"),
               stringsAsFactors = FALSE)

month_levels <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

rain <- df |>
  select(YEAR, all_of(month_levels)) |>
  pivot_longer(-YEAR, names_to = "month_abbr", values_to = "precip") |>
  mutate(
    precip    = na_if(precip, -999),
    month_num = match(.data$month_abbr, month_levels),
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

print(rain |> features(precip, feat_stl))

# Statistical tests
cat("\n== ADF test ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS test ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box test: raw series ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

cat("\n== Mann-Kendall trend test ==\n")
print(Kendall::MannKendall(rain$precip))

cat("\n== Minimum precipitation ==\n")
print(min(rain$precip, na.rm = TRUE))

# Train/test split
h     <- 12
train <- rain |> filter(month <= max(month) - h)

# Model
fit <- train |> model(
  snaive      = SNAIVE(precip),
  arima_four4 = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0))
)
fc <- fit |> forecast(h = h)

# Differenced-series diagnostics
arima_d <- train |>
  features(precip, unitroot_ndiffs) |>
  pull(ndiffs)

arima_series <- if (arima_d > 0) {
  train |> mutate(precip_diff = difference(precip, differences = arima_d))
} else {
  train |> mutate(precip_diff = precip)
}

p_acf_diff <- arima_series |>
  ACF(precip_diff, lag_max = 36) |>
  autoplot() +
  labs(title = paste0("ACF after differencing (d = ", arima_d, ")"))

p_pacf_diff <- arima_series |>
  PACF(precip_diff, lag_max = 36) |>
  autoplot() +
  labs(title = paste0("PACF after differencing (d = ", arima_d, ")"))

p_arima_correlation <- patchwork::wrap_plots(p_acf_diff, p_pacf_diff, ncol = 2)

# Accuracy
accuracy_display <- bind_rows(
  fit |> accuracy() |> mutate(data_set = "Training set"),
  fc |> accuracy(rain) |> mutate(data_set = "Test set")
) |>
  transmute(
    model = recode(.model,
                   snaive = "Seasonal naive",
                   arima_four4 = "ARIMA + Fourier"),
    data_set, RMSE, MAE, MAPE, MASE
  ) |>
  arrange(model, factor(data_set, levels = c("Training set", "Test set"))) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

cat("\n--- ARIMA and Benchmark Accuracy ---\n")
print(accuracy_display)
fit |> select(arima_four4) |> report()
fit |> select(arima_four4) |> gg_tsresiduals()

# Residual checks
arima_coefs  <- fit |> select(arima_four4) |> tidy()
arima_dof    <- sum(grepl("^(ar|ma)[0-9]+$", arima_coefs$term))

cat("\n== Ljung-Box test: ARIMA residuals ==\n")
print(augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 12, dof = arima_dof))
print(augment(fit) |> filter(.model == "arima_four4") |> features(.innov, ljung_box, lag = 24, dof = arima_dof))

cat("\n== Residual ACF summary ==\n")
print(augment(fit) |> filter(.model == "arima_four4") |> as_tibble() |>
        summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
                  n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)))

# Train/test comparison
acc_train <- fit |> accuracy() |> filter(.model == "arima_four4") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test  <- fc |> accuracy(rain) |> filter(.model == "arima_four4") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
holdout <- acc_train |> left_join(acc_test, by = ".model") |>
  select(.model, MASE_train, RMSE_train, MASE_test, RMSE_test)

# Rolling-origin CV
n_folds_clean <- length(seq(360, nrow(train) - h, by = 6))
cv_data <- train |>
  stretch_tsibble(.init = 360, .step = 6) |>
  filter(.id <= n_folds_clean)

set.seed(2026)
cv_acc <- cv_data |>
  model(arima_four4 = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0))) |>
  forecast(h = h) |>
  accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  model         = "arima_fourier",
  MASE_train    = holdout$MASE_train,
  RMSE_train    = holdout$RMSE_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MASE_cv    = cv_summary$sd_MASE,
  n_folds       = cv_summary$n_folds
)

# Output
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/arima_results.csv", row.names = FALSE)
write.csv(accuracy_display, "output/arima_accuracy_train_test.csv", row.names = FALSE)
ggsave("output/plots/group_summary/arima_differenced_acf_pacf.png",
       p_arima_correlation, width = 10, height = 4.5, dpi = 150)

# Forecast plot
plot_from       <- yearmonth("2013 Jan")
test_actual_tbl <- rain |> as_tibble() |> filter(month > max(train$month)) |>
  transmute(month, precip)

p_fc <- fc |> filter(.model == "arima_four4") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  geom_line(data = test_actual_tbl, aes(x = month, y = precip),
            color = "red", linewidth = 0.45) +
  labs(title = "ARIMA + Fourier(K=4): Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
ggsave("output/plots/group_summary/fc_arima.png", p_fc, width = 8, height = 5, dpi = 150)

ggsave("output/plots/group_summary/resid_arima.png",
       fit |> select(arima_four4) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)

cat("\nDone. Wrote ARIMA result and train/test accuracy CSV files, plus 3 plots.\n")
