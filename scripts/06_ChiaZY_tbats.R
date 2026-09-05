# TBATS rainfall forecast

# Setup
pkgs <- c("fpp3", "tseries", "zoo", "forecast", "Kendall")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
library(purrr)
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

# STL decomposition
p_stl <- rain |> model(STL(precip)) |> components() |> autoplot() +
  labs(title = "STL decomposition")
print(p_stl)
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
h        <- 12
train    <- rain |> filter(month <= max(month) - h)
train_ts <- ts(train$precip, frequency = 12)

# Model
fit_tbats <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                             use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats  <- forecast::forecast(fit_tbats, h = h)
print(fit_tbats)
cat("\nBox-Cox lambda:", if (is.null(fit_tbats$lambda)) "none" else fit_tbats$lambda, "\n")

test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats   <- forecast::accuracy(
  fc_tbats,
  test_actual,
  d = 0,
  D = 1
)

accuracy_display <- acc_tbats |>
  as.data.frame() |>
  tibble::rownames_to_column("data_set") |>
  as_tibble() |>
  transmute(model = "TBATS", data_set, RMSE, MAE, MAPE, MASE) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

snaive_fit <- train |> model(snaive = SNAIVE(precip))
snaive_fc  <- snaive_fit |> forecast(h = h)
snaive_accuracy <- bind_rows(
  snaive_fit |> accuracy() |> mutate(data_set = "Training set"),
  snaive_fc |> accuracy(rain) |> mutate(data_set = "Test set")
) |>
  transmute(model = "Seasonal naive", data_set, RMSE, MAE, MAPE, MASE) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

accuracy_display <- bind_rows(accuracy_display, snaive_accuracy) |>
  arrange(model, factor(data_set, levels = c("Training set", "Test set")))

cat("\n--- TBATS and Benchmark Accuracy ---\n")
print(accuracy_display)

# Residual checks
resid_tbats <- residuals(fit_tbats)

cat("\n== Ljung-Box test: TBATS residuals ==\n")
print(Box.test(resid_tbats, lag = 12, type = "Ljung-Box"))
print(Box.test(resid_tbats, lag = 24, type = "Ljung-Box"))

cat("\n== Residual ACF summary ==\n")
cat("n_lags_out_12:", acf_out_of_bounds(resid_tbats, lag.max = 12),
    " n_lags_out_24:", acf_out_of_bounds(resid_tbats, lag.max = 24), "\n")

# Train/test comparison
mase_train <- acc_tbats["Training set", "MASE"]
rmse_train <- acc_tbats["Training set", "RMSE"]

# Rolling-origin CV
origins <- seq(360, nrow(train) - h, by = 6)

cv_tbats <- map_dfr(origins, function(i) {
  tr  <- train |> slice(1:i)
  te  <- train |> slice((i + 1):(i + h)) |> pull(precip)
  m   <- forecast::tbats(ts(tr$precip, frequency = 12), use.box.cox = NULL,
                         use.trend = FALSE, use.damped.trend = FALSE,
                         seasonal.periods = 12)
  fcv <- forecast::forecast(m, h = h)
  acc <- forecast::accuracy(
    fcv,
    te,
    d = 0,
    D = 1
  )
  tibble(MASE = acc["Test set", "MASE"], RMSE = acc["Test set", "RMSE"])
})

cv_summary <- cv_tbats |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  model         = "tbats",
  MASE_train    = mase_train,
  RMSE_train    = rmse_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MASE_cv    = cv_summary$sd_MASE,
  n_folds       = cv_summary$n_folds,
  lb_pvalue_24  = Box.test(resid_tbats, lag = 24, type = "Ljung-Box")$p.value
)

# Output
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/tbats_results.csv", row.names = FALSE)
write.csv(accuracy_display, "output/tbats_accuracy_train_test.csv", row.names = FALSE)

# STL plot
ggsave("output/plots/group_summary/stl_tbats_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)

# Forecast plot
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
ggsave("output/plots/group_summary/fc_tbats.png", p_fc, width = 8, height = 5, dpi = 150)

png("output/plots/group_summary/resid_tbats.png", width = 800, height = 600, res = 150)
forecast::checkresiduals(fit_tbats)
dev.off()

cat("\nDone. Wrote TBATS result and train/test accuracy CSV files, plus 3 plots.\n")
