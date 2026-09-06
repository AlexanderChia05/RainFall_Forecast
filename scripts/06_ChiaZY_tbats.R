# TBATS rainfall forecast

# Shared setup and data acquisition make this script independently runnable.
source("scripts/00_setup.R")
source("scripts/01_data_pull.R")

cat("rain:", nrow(rain), "obs,", format(min(rain$month)), "to", format(max(rain$month)), "\n")
cat("Missing before interpolation:", missing_before,
    "| Missing after interpolation:", missing_after, "\n")

# EDA
p_raw <- rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
p_season <- rain |> gg_season(precip) + labs(title = "Seasonal plot - monsoon cycle")
p_subseries <- rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
p_acf <- rain |> ACF(precip, lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
p_pacf <- rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# STL decomposition
p_stl <- rain |> model(STL(precip)) |> components() |> autoplot() +
  labs(title = "STL decomposition")
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

tbats_specification <- tibble(
  model = "TBATS",
  box_cox_selection = "automatic by AIC",
  lambda = if (is.null(fit_tbats$lambda)) NA_real_ else fit_tbats$lambda,
  ar_order = length(fit_tbats$ar.coefficients),
  ma_order = length(fit_tbats$ma.coefficients),
  trend = !is.null(fit_tbats$beta),
  damped_trend = !is.null(fit_tbats$damping.parameter),
  seasonal_period = paste(fit_tbats$seasonal.periods, collapse = ","),
  harmonics = paste(fit_tbats$k.vector, collapse = ","),
  mase_d = 0L,
  mase_D = 1L,
  ljung_box_fitdf = 0L
)

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
dir.create("output/tables/model_details", recursive = TRUE, showWarnings = FALSE)
write.csv(tbats_specification, "output/tables/model_details/tbats_parameters.csv", row.names = FALSE)

# STL plot
suppressMessages(
  suppressWarnings(
    ggsave("output/plots/group_summary/stl_tbats_decomposition.png", p_stl,
           width = 9, height = 6, dpi = 150)
  )
)

# Standardized forecast and residual plots shared with the other models.
save_tbats_forecast_plot(
  fc_tbats, rain, train, "output/plots/group_summary/fc_tbats.png"
)
save_residual_diagnostic(
  resid_tbats, train$month, "TBATS",
  "output/plots/group_summary/resid_tbats.png"
)

cat("\nDone. Wrote the TBATS parameter table and 3 plots.\n")
