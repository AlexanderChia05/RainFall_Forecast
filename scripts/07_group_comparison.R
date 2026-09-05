source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

h <- 12
train <- rain |> filter(month <= max(month) - h)

# Models
set.seed(2026)
fit <- train |> model(
  seasonal_naive = SNAIVE(precip),
  ets_additive   = ETS(precip ~ error("A") + trend("N") + season("A")),
  arima_fourier  = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
  tslm_fourier   = TSLM(precip ~ trend() + fourier(K = 4))
)

fc <- fit |> forecast(h = h)
acc_test <- fc |> accuracy(rain) |>
  select(model = .model, MASE, RMSE, MAE, MAPE)

# Residual checks
arima_coefs <- fit |> select(arima_fourier) |> tidy()
arima_dof   <- sum(grepl("^(ar|ma)[0-9]+$", arima_coefs$term))

lb12 <- bind_rows(
  augment(fit) |> filter(.model != "arima_fourier") |> features(.innov, ljung_box, lag = 12),
  augment(fit) |> filter(.model == "arima_fourier") |> features(.innov, ljung_box, lag = 12, dof = arima_dof)
) |> rename(model = .model, lb_stat_12 = lb_stat, lb_pvalue_12 = lb_pvalue)

lb24 <- bind_rows(
  augment(fit) |> filter(.model != "arima_fourier") |> features(.innov, ljung_box, lag = 24),
  augment(fit) |> filter(.model == "arima_fourier") |> features(.innov, ljung_box, lag = 24, dof = arima_dof)
) |> rename(model = .model, lb_stat_24 = lb_stat, lb_pvalue_24 = lb_pvalue)

acf_check <- fit |>
  augment() |>
  as_tibble() |>
  group_by(.model) |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)) |>
  rename(model = .model)

acc_train <- fit |> accuracy() |>
  select(model = .model,
         MASE_train = MASE, RMSE_train = RMSE,
         MAE_train = MAE, MAPE_train = MAPE)

# TBATS
train_ts   <- ts(train$precip, frequency = 12)
fit_tbats  <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                               use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats   <- forecast::forecast(fit_tbats, h = h)
test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats_full <- forecast::accuracy(
  fc_tbats,
  test_actual,
  d = 0,
  D = 1
)
resid_tbats <- residuals(fit_tbats)

lb_tbats12  <- Box.test(resid_tbats, lag = 12, type = "Ljung-Box")
lb_tbats24  <- Box.test(resid_tbats, lag = 24, type = "Ljung-Box")

acc_test <- acc_test |> bind_rows(tibble(
  model = "tbats",
  MASE = acc_tbats_full["Test set", "MASE"], RMSE = acc_tbats_full["Test set", "RMSE"],
  MAE  = acc_tbats_full["Test set", "MAE"],  MAPE = acc_tbats_full["Test set", "MAPE"]
)) |> arrange(MASE)

lb12 <- lb12 |> bind_rows(tibble(
  model = "tbats", lb_stat_12 = unname(lb_tbats12$statistic), lb_pvalue_12 = lb_tbats12$p.value
))
lb24 <- lb24 |> bind_rows(tibble(
  model = "tbats", lb_stat_24 = unname(lb_tbats24$statistic), lb_pvalue_24 = lb_tbats24$p.value
))

acf_check <- acf_check |> bind_rows(tibble(
  model = "tbats",
  n_lags_out_12 = acf_out_of_bounds(resid_tbats, lag.max = 12),
  n_lags_out_24 = acf_out_of_bounds(resid_tbats, lag.max = 24)
))

acc_train <- acc_train |> bind_rows(tibble(
  model = "tbats",
  MASE_train = acc_tbats_full["Training set", "MASE"],
  RMSE_train = acc_tbats_full["Training set", "RMSE"],
  MAE_train  = acc_tbats_full["Training set", "MAE"],
  MAPE_train = acc_tbats_full["Training set", "MAPE"]
))

# Rolling-origin CV
n_folds_clean <- length(seq(360, nrow(train) - h, by = 6))
origins <- seq(360, nrow(train) - h, by = 6)
cv_data <- train |>
  stretch_tsibble(.init = 360, .step = 6) |>
  filter(.id <= n_folds_clean)

set.seed(2026)
cv_fits <- cv_data |>
  model(
    seasonal_naive = SNAIVE(precip),
    ets_additive   = ETS(precip ~ error("A") + trend("N") + season("A")),
    arima_fourier  = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
    tslm_fourier   = TSLM(precip ~ trend() + fourier(K = 4))
  )

cv_acc <- cv_fits |> forecast(h = h) |> accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  group_by(model = .model) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE = min(MASE), max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds = n())

# TBATS CV
cv_tbats <- map_dfr(origins, function(i) {
  tr  <- train |> slice(1:i)
  te  <- train |> slice((i + 1):(i + h)) |> pull(precip)
  m   <- forecast::tbats(ts(tr$precip, frequency = 12), use.box.cox = NULL,
                          use.trend = FALSE, use.damped.trend = FALSE,
                          seasonal.periods = 12)
  fc  <- forecast::forecast(m, h = h)
  acc <- forecast::accuracy(
    fc,
    te,
    d = 0,
    D = 1
  )
  tibble(MASE = acc["Test set", "MASE"], RMSE = acc["Test set", "RMSE"])
})
cv_summary <- cv_summary |> bind_rows(
  cv_tbats |> summarise(model = "tbats", mean_MASE = mean(MASE), sd_MASE = sd(MASE),
                         min_MASE = min(MASE), max_MASE = max(MASE),
                         mean_RMSE = mean(RMSE), n_folds = n())
) |> arrange(mean_MASE)
write.csv(cv_summary, "output/model_comparison_cv_summary.csv", row.names = FALSE)

summary_tbl <- acc_test |>
  left_join(lb12, by = "model") |>
  left_join(lb24, by = "model") |>
  left_join(acf_check, by = "model")

# Display tables
accuracy_comparison <- bind_rows(
  acc_train |>
    transmute(model, data_set = "Training set",
              RMSE = RMSE_train, MAE = MAE_train,
              MAPE = MAPE_train, MASE = MASE_train),
  acc_test |>
    transmute(model, data_set = "Test set", RMSE, MAE, MAPE, MASE)
) |>
  mutate(model = recode(
    model,
    seasonal_naive = "Seasonal naive",
    ets_additive = "ETS",
    arima_fourier = "ARIMA + Fourier",
    tslm_fourier = "TSLM + Fourier",
    tbats = "TBATS"
  )) |>
  arrange(model, factor(data_set, levels = c("Training set", "Test set"))) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

cv_comparison <- cv_summary |>
  transmute(
    model = recode(
      model,
      seasonal_naive = "Seasonal naive",
      ets_additive = "ETS",
      arima_fourier = "ARIMA + Fourier",
      tslm_fourier = "TSLM + Fourier",
      tbats = "TBATS"
    ),
    CV_RMSE = round(mean_RMSE, 3),
    CV_MASE = round(mean_MASE, 3),
    SD_MASE = round(sd_MASE, 3),
    folds = n_folds
  ) |>
  arrange(CV_MASE)

diagnostics_comparison <- lb12 |>
  left_join(lb24, by = "model") |>
  left_join(acf_check, by = "model") |>
  transmute(
    model = recode(
      model,
      seasonal_naive = "Seasonal naive",
      ets_additive = "ETS",
      arima_fourier = "ARIMA + Fourier",
      tslm_fourier = "TSLM + Fourier",
      tbats = "TBATS"
    ),
    LB_pvalue_12 = round(lb_pvalue_12, 3),
    LB_pvalue_24 = round(lb_pvalue_24, 3),
    ACF_lags_out_12 = n_lags_out_12,
    ACF_lags_out_24 = n_lags_out_24
  )

cat("\n--- Training and Test Accuracy ---\n")
print(accuracy_comparison)
cat("\n--- Rolling-origin Cross-validation ---\n")
print(cv_comparison)
cat("\n--- Residual Diagnostics ---\n")
print(diagnostics_comparison)

dir.create("output", showWarnings = FALSE)
write.csv(summary_tbl, "output/model_comparison_summary.csv", row.names = FALSE)
write.csv(accuracy_comparison, "output/model_accuracy_train_test.csv", row.names = FALSE)
write.csv(cv_comparison, "output/model_cv_comparison.csv", row.names = FALSE)
write.csv(diagnostics_comparison, "output/model_diagnostics.csv", row.names = FALSE)
