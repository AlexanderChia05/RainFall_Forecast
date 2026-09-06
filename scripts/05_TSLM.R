# TSLM with trend and Fourier terms rainfall forecast

# Shared setup and data acquisition make this script independently runnable.
source("scripts/00_setup.R")
if (!file.exists("data/rain.rds")) source("scripts/01_data_pull.R")
rain <- readRDS("data/rain.rds")

cat("rain:", nrow(rain), "obs,", format(min(rain$month)), "to", format(max(rain$month)),
    "| missing:", sum(is.na(rain$precip)), "\n")

# Dataset characteristics quoted in the individual report. They are repeated here so
# that this script reproduces every number it cites without sourcing another model
# script or 02_eda_stationarity.R.
cat("\n== STL feature strengths ==\n")
print(rain |> features(precip, feat_stl))

cat("\n== ADF test ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS test ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box test: raw series ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

cat("\n== Mann-Kendall trend test ==\n")
print(Kendall::MannKendall(rain$precip))

cat("\n== Distribution of the raw series ==\n")
print(rain |> as_tibble() |> summarise(
  n    = n(),          missing = sum(is.na(precip)),
  min  = min(precip),  max     = max(precip),
  mean = mean(precip), sd      = sd(precip)
))

# Train/test split
h     <- 12
train <- rain |> filter(month <= max(month) - h)

# Model
fit <- train |> model(
  snaive = SNAIVE(precip),
  tslm   = TSLM(precip ~ trend() + fourier(K = 3))
)
fc <- fit |> forecast(h = h)

# Accuracy
accuracy_display <- bind_rows(
  fit |> accuracy() |> mutate(data_set = "Training set"),
  fc |> accuracy(rain) |> mutate(data_set = "Test set")
) |>
  transmute(
    model = recode(.model, snaive = "Seasonal naive", tslm = "TSLM + Fourier"),
    data_set, RMSE, MAE, MAPE, MASE
  ) |>
  arrange(model, factor(data_set, levels = c("Training set", "Test set"))) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

cat("\n--- TSLM and Benchmark Accuracy ---\n")
print(accuracy_display)

fit |> select(tslm) |> report()
fit |> select(tslm) |> gg_tsresiduals()

tslm_parameters <- fit |>
  select(tslm) |>
  tidy() |>
  mutate(model = "TSLM with trend and Fourier K=3") |>
  relocate(model)

tslm_specification <- tibble(
  model = "TSLM with trend and Fourier terms",
  estimation = "ordinary least squares",
  trend = "linear",
  seasonal_period = 12L,
  fourier_K = 3L,
  fourier_regressors = 8L,
  autoregressive_errors = FALSE
)

# Residual checks
cat("\n== Ljung-Box test: TSLM residuals ==\n")
print(augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 12))
print(augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 24))

cat("\n== Residual ACF summary ==\n")
print(augment(fit) |> filter(.model == "tslm") |> as_tibble() |>
        summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
                  n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)))

# Train/test comparison
acc_train <- fit |> accuracy() |> filter(.model == "tslm") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test  <- fc |> accuracy(rain) |> filter(.model == "tslm") |>
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
  model(tslm = TSLM(precip ~ trend() + fourier(K = 3))) |>
  forecast(h = h) |>
  accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  model         = "tslm_fourier",
  MASE_train    = holdout$MASE_train,
  RMSE_train    = holdout$RMSE_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MASE_cv    = cv_summary$sd_MASE,
  n_folds       = cv_summary$n_folds
)

# Output
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables/model_details", recursive = TRUE, showWarnings = FALSE)
write.csv(tslm_parameters, "output/tables/model_details/tslm_parameters.csv", row.names = FALSE)

# Forecast and residual plots
save_fable_forecast_plot(
  fc |> filter(.model == "tslm"), rain, "TSLM + Fourier",
  "TSLM (trend + Fourier K=3): Forecast vs Actual",
  "output/plots/group_summary/fc_tslm.png"
)
u_tslm <- augment(fit) |> filter(.model == "tslm") |> as_tibble()
save_residual_diagnostic(u_tslm$.innov, u_tslm$month, "TSLM + Fourier",
                         "output/plots/group_summary/resid_tslm.png")

cat("\nDone. Wrote the TSLM parameter table and 2 plots.\n")
