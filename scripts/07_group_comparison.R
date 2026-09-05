source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

h <- 12
train <- rain |> filter(month <= max(month) - h)

# The 3 fable-native models
set.seed(2026)
fit <- train |> model(
  snaive  = SNAIVE(precip),
  ets_A   = ETS(precip ~ error("A") + trend("N") + season("A")),
  arima_B = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
  tslm_C  = TSLM(precip ~ trend() + fourier(K = 4))
)

fc <- fit |> forecast(h = h)
acc_test <- fc |> accuracy(rain) |>
  select(member = .model, MASE, RMSE, MAE, MAPE)

# ARIMA's Ljung-Box degrees of freedom must be reduced by the number of
# AR/MA parameters it estimated, or the test overstates residual
# randomness. Counted from the fitted coefficient names rather than
# assumed - PDQ(0,0,0) means no seasonal ar/ma term to also count, and
# the Fourier regressor coefficients are named after the regressor, not
# ar/ma, so this pattern only ever matches genuine ARMA error terms.
# ETS and TSLM keep the default dof = 0 as basic residual screening.
arima_coefs <- fit |> select(arima_B) |> tidy()
arima_dof   <- sum(grepl("^(ar|ma)[0-9]+$", arima_coefs$term))
cat("ARIMA AR/MA parameter count used as Ljung-Box dof:", arima_dof, "\n")

lb12 <- bind_rows(
  augment(fit) |> filter(.model != "arima_B") |> features(.innov, ljung_box, lag = 12),
  augment(fit) |> filter(.model == "arima_B") |> features(.innov, ljung_box, lag = 12, dof = arima_dof)
) |> rename(member = .model, lb_stat_12 = lb_stat, lb_pvalue_12 = lb_pvalue)

lb24 <- bind_rows(
  augment(fit) |> filter(.model != "arima_B") |> features(.innov, ljung_box, lag = 24),
  augment(fit) |> filter(.model == "arima_B") |> features(.innov, ljung_box, lag = 24, dof = arima_dof)
) |> rename(member = .model, lb_stat_24 = lb_stat, lb_pvalue_24 = lb_pvalue)

acf_check <- fit |>
  augment() |>
  as_tibble() |>
  group_by(.model) |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)) |>
  rename(member = .model)

acc_train <- fit |> accuracy() |>
  select(member = .model, MASE_train = MASE, RMSE_train = RMSE)

# TBATS
train_ts   <- ts(train$precip, frequency = 12)
fit_tbats  <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                               use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats   <- forecast::forecast(fit_tbats, h = h)
test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats_full <- forecast::accuracy(fc_tbats, test_actual)
resid_tbats <- residuals(fit_tbats)

# TBATS's Ljung-Box fitdf must likewise account for the parameters it
# estimated internally (Box-Cox, ARMA errors, seasonal states);
# Box.test()'s fitdf defaults to 0, so it is set explicitly here.
tbats_dof   <- forecast::modeldf(fit_tbats)
cat("TBATS model degrees of freedom used as Ljung-Box fitdf:", tbats_dof, "\n")
lb_tbats12  <- Box.test(resid_tbats, lag = 12, type = "Ljung-Box", fitdf = tbats_dof)
lb_tbats24  <- Box.test(resid_tbats, lag = 24, type = "Ljung-Box", fitdf = tbats_dof)

acc_test <- acc_test |> bind_rows(tibble(
  member = "tbats_D",
  MASE = acc_tbats_full["Test set", "MASE"], RMSE = acc_tbats_full["Test set", "RMSE"],
  MAE  = acc_tbats_full["Test set", "MAE"],  MAPE = acc_tbats_full["Test set", "MAPE"]
)) |> arrange(MASE)
print(acc_test)

lb12 <- lb12 |> bind_rows(tibble(
  member = "tbats_D", lb_stat_12 = unname(lb_tbats12$statistic), lb_pvalue_12 = lb_tbats12$p.value
))
lb24 <- lb24 |> bind_rows(tibble(
  member = "tbats_D", lb_stat_24 = unname(lb_tbats24$statistic), lb_pvalue_24 = lb_tbats24$p.value
))
print(lb12)
print(lb24)

acf_check <- acf_check |> bind_rows(tibble(
  member = "tbats_D",
  n_lags_out_12 = acf_out_of_bounds(resid_tbats, lag.max = 12),
  n_lags_out_24 = acf_out_of_bounds(resid_tbats, lag.max = 24)
))
print(acf_check)

acc_train <- acc_train |> bind_rows(tibble(
  member = "tbats_D",
  MASE_train = acc_tbats_full["Training set", "MASE"],
  RMSE_train = acc_tbats_full["Training set", "RMSE"]
))

# Rolling-origin CV (.init=360, .step=6). Folds run only over the
# training window (up to 2024-12), never the 2025 holdout - otherwise
# the same 2025 observations that back the holdout accuracy above would
# also leak into the CV folds and inflate both numbers on the same data.
n_folds_clean <- length(seq(360, nrow(train) - h, by = 6))
origins <- seq(360, nrow(train) - h, by = 6)

set.seed(2026)
cv_fits <- train |>
  stretch_tsibble(.init = 360, .step = 6) |>
  model(
    snaive  = SNAIVE(precip),
    ets_A   = ETS(precip ~ error("A") + trend("N") + season("A")),
    arima_B = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
    tslm_C  = TSLM(precip ~ trend() + fourier(K = 4))
  )

cv_acc <- cv_fits |> forecast(h = 12) |> accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  filter(!is.na(MASE), .id <= n_folds_clean) |>
  group_by(member = .model) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE = min(MASE), max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds = n())

# TBATS CV - every slice drawn from train, not rain, so 2025 is never
# touched inside this loop either.
cv_tbats <- map_dfr(origins, function(i) {
  tr  <- train |> slice(1:i)
  te  <- train |> slice((i + 1):(i + h)) |> pull(precip)
  m   <- forecast::tbats(ts(tr$precip, frequency = 12), use.box.cox = NULL,
                          use.trend = FALSE, use.damped.trend = FALSE,
                          seasonal.periods = 12)
  fc  <- forecast::forecast(m, h = h)
  acc <- forecast::accuracy(fc, te)
  tibble(MASE = acc["Test set", "MASE"], RMSE = acc["Test set", "RMSE"])
})
cv_summary <- cv_summary |> bind_rows(
  cv_tbats |> summarise(member = "tbats_D", mean_MASE = mean(MASE), sd_MASE = sd(MASE),
                         min_MASE = min(MASE), max_MASE = max(MASE),
                         mean_RMSE = mean(RMSE), n_folds = n())
) |> arrange(mean_MASE)
print(cv_summary)
write.csv(cv_summary, "output/model_comparison_cv_summary.csv", row.names = FALSE)

# Overfitting check
gap_tbl <- acc_train |>
  left_join(acc_test |> select(member, MASE_test_holdout = MASE, RMSE_test_holdout = RMSE), by = "member") |>
  left_join(cv_summary |> select(member, MASE_test_cv = mean_MASE, RMSE_test_cv = mean_RMSE), by = "member") |>
  mutate(gap_pct_holdout = abs(MASE_test_holdout - MASE_train) / MASE_test_holdout,
         gap_pct_cv       = abs(MASE_test_cv - MASE_train) / MASE_test_cv,
         rmse_ratio_holdout = RMSE_test_holdout / RMSE_train,
         rmse_ratio_cv       = RMSE_test_cv / RMSE_train,
         within_10pct_cv  = gap_pct_cv <= 0.10,
         within_1_3x_cv   = rmse_ratio_cv <= 1.3) |>
  arrange(rmse_ratio_cv)
print(gap_tbl)

summary_tbl <- acc_test |>
  left_join(lb12, by = "member") |>
  left_join(lb24, by = "member") |>
  left_join(acf_check, by = "member") |>
  left_join(gap_tbl |> select(member, gap_pct_holdout, gap_pct_cv, rmse_ratio_holdout,
                               rmse_ratio_cv, within_10pct_cv, within_1_3x_cv), by = "member")
print(summary_tbl)

dir.create("output", showWarnings = FALSE)
write.csv(summary_tbl, "output/model_comparison_summary.csv", row.names = FALSE)
