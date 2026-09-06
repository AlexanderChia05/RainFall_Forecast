source("scripts/00_setup.R")
if (!file.exists("data/rain.rds")) source("scripts/01_data_pull.R")
rain <- readRDS("data/rain.rds")

h <- 12
train <- rain |> filter(month <= max(month) - h)

# Models
set.seed(2026)
fit <- train |> model(
  seasonal_naive = SNAIVE(precip),
  ets_additive   = ETS(precip ~ error("A") + trend("N") + season("A")),
  arima_fourier  = ARIMA(precip ~ fourier(K = 3) + pdq() + PDQ(0, 0, 0)),
  tslm_fourier   = TSLM(precip ~ trend() + fourier(K = 3))
)

fc <- fit |> forecast(h = h)
acc_test <- fc |> accuracy(rain) |>
  select(model = .model, MASE, RMSE, MAE, MAPE)

# Residual checks
arima_coefs <- fit |> select(arima_fourier) |> tidy()
arima_dof   <- sum(grepl("^(ar|ma)[0-9]+$", arima_coefs$term))
arima_p     <- sum(grepl("^ar[0-9]+$", arima_coefs$term))
arima_q     <- sum(grepl("^ma[0-9]+$", arima_coefs$term))
arima_d     <- train |>
  features(precip, unitroot_ndiffs) |>
  pull(ndiffs)

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
    arima_fourier  = ARIMA(precip ~ fourier(K = 3) + pdq() + PDQ(0, 0, 0)),
    tslm_fourier   = TSLM(precip ~ trend() + fourier(K = 3))
  )

cv_acc <- cv_fits |> forecast(h = h) |> accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  group_by(model = .model) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE = min(MASE), max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), mean_MAE = mean(MAE),
            mean_MAPE = mean(MAPE), n_folds = n())

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
  tibble(
    MASE = acc["Test set", "MASE"], RMSE = acc["Test set", "RMSE"],
    MAE = acc["Test set", "MAE"], MAPE = acc["Test set", "MAPE"],
    coverage_80 = mean(te >= fc$lower[, "80%"] & te <= fc$upper[, "80%"]),
    coverage_95 = mean(te >= fc$lower[, "95%"] & te <= fc$upper[, "95%"]),
    width_80 = mean(fc$upper[, "80%"] - fc$lower[, "80%"]),
    width_95 = mean(fc$upper[, "95%"] - fc$lower[, "95%"])
  )
})

tbats_interval_calibration <- tibble(
  nominal_level = c(80L, 95L),
  empirical_coverage = c(mean(cv_tbats$coverage_80), mean(cv_tbats$coverage_95)),
  mean_interval_width = c(mean(cv_tbats$width_80), mean(cv_tbats$width_95)),
  forecast_count = nrow(cv_tbats) * h
)
cv_summary <- cv_summary |> bind_rows(
  cv_tbats |> summarise(model = "tbats", mean_MASE = mean(MASE), sd_MASE = sd(MASE),
                         min_MASE = min(MASE), max_MASE = max(MASE),
                         mean_RMSE = mean(RMSE), mean_MAE = mean(MAE),
                         mean_MAPE = mean(MAPE), n_folds = n())
) |> arrange(mean_MASE)

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
    CV_MAE = round(mean_MAE, 3),
    CV_MAPE = round(mean_MAPE, 3),
    CV_MASE = round(mean_MASE, 3),
    SD_MASE = round(sd_MASE, 3),
    Min_MASE = round(min_MASE, 3),
    Max_MASE = round(max_MASE, 3),
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

tbats_model_label <- paste0(
  "TBATS(",
  if (is.null(fit_tbats$lambda)) "-" else round(fit_tbats$lambda, 3),
  ", {", length(fit_tbats$ar.coefficients), ",", length(fit_tbats$ma.coefficients), "}, ",
  if (is.null(fit_tbats$damping.parameter)) "-" else round(fit_tbats$damping.parameter, 3),
  ", {<", paste(fit_tbats$seasonal.periods, collapse = ","), ",",
  paste(fit_tbats$k.vector, collapse = ","), ">})"
)

model_specifications <- tibble(
  model = c("Seasonal naive", "ETS", "ARIMA + Fourier", "TSLM + Fourier", "TBATS"),
  specification = c(
    "SNAIVE with seasonal period 12",
    "ETS(A,N,A)",
    paste0("ARIMA(", arima_p, ",", arima_d, ",", arima_q,
           ") errors; seasonal order (0,0,0); Fourier K=3"),
    "Linear trend; Fourier K=3; ordinary least squares",
    tbats_model_label
  ),
  seasonal_period = 12L,
  mase_scaling = "seasonal naive at lag 12",
  ljung_box_fitdf = c(0L, 0L, arima_dof, 0L, 0L)
)

package_versions <- tibble(
  component = c("R", "fpp3", "fable", "forecast", "feasts", "tsibble"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(c("fpp3", "fable", "forecast", "feasts", "tsibble"),
           function(pkg) as.character(packageVersion(pkg)), character(1))
  )
)

model_results <- accuracy_comparison |>
  mutate(data_set = recode(data_set,
                           "Training set" = "Training",
                           "Test set" = "Test")) |>
  pivot_wider(
    names_from = data_set,
    values_from = c(RMSE, MAE, MAPE, MASE),
    names_glue = "{data_set}_{.value}"
  ) |>
  left_join(cv_comparison, by = "model") |>
  left_join(diagnostics_comparison, by = "model") |>
  select(
    model,
    Training_RMSE, Training_MAE, Training_MAPE, Training_MASE,
    Test_RMSE, Test_MAE, Test_MAPE, Test_MASE,
    CV_RMSE, CV_MAE, CV_MAPE, CV_MASE, SD_MASE, Min_MASE, Max_MASE, folds,
    LB_pvalue_12, LB_pvalue_24, ACF_lags_out_12, ACF_lags_out_24
  )

cat("\n--- Training and Test Accuracy ---\n")
print(accuracy_comparison)
cat("\n--- Rolling-origin Cross-validation ---\n")
print(cv_comparison)
cat("\n--- Residual Diagnostics ---\n")
print(diagnostics_comparison)
cat("\n--- Model Specifications ---\n")
print(model_specifications)
cat("\n--- TBATS Interval Calibration ---\n")
print(tbats_interval_calibration)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(model_results, "output/tables/model_results.csv", row.names = FALSE)
write.csv(model_specifications, "output/tables/model_specifications.csv", row.names = FALSE)
write.csv(package_versions, "output/tables/package_versions.csv", row.names = FALSE)
write.csv(
  tbats_interval_calibration,
  "output/tables/tbats_interval_calibration.csv",
  row.names = FALSE
)

# Candidate-specification audit --------------------------------------------
# Select each model family's specification using the same rolling origins,
# then use the untouched 2025 holdout as a separate confirmation.
run_model_selection_audit <- function() {
score_fable_cv <- function(candidate_fit, label) {
  candidate_fit |>
    forecast(h = h) |>
    accuracy(train, by = c(".model", ".id")) |>
    summarise(
      specification = label,
      CV_MASE = mean(MASE),
      CV_RMSE = mean(RMSE),
      CV_MAE = mean(MAE),
      CV_MAPE = mean(MAPE),
      SD_MASE = sd(MASE),
      Min_MASE = min(MASE),
      Max_MASE = max(MASE),
      folds = n(),
      .groups = "drop"
    )
}

score_fable_holdout <- function(candidate_fit, label) {
  candidate_fit |>
    forecast(h = h) |>
    accuracy(rain) |>
    transmute(
      specification = label,
      Test_MASE = MASE,
      Test_RMSE = RMSE,
      Test_MAE = MAE,
      Test_MAPE = MAPE
    )
}

ets_audit <- bind_rows(
  left_join(
    score_fable_cv(
      cv_data |> model(candidate = ETS(precip ~ error("A") + trend("N") + season("A"))),
      "ETS(A,N,A)"
    ),
    score_fable_holdout(
      train |> model(candidate = ETS(precip ~ error("A") + trend("N") + season("A"))),
      "ETS(A,N,A)"
    ),
    by = "specification"
  ),
  left_join(
    score_fable_cv(cv_data |> model(candidate = ETS(precip)), "ETS automatic"),
    score_fable_holdout(train |> model(candidate = ETS(precip)), "ETS automatic"),
    by = "specification"
  )
)

fourier_audit <- purrr::map_dfr(1:5, function(k) {
  arima_label <- paste0("ARIMA + Fourier K=", k)
  tslm_label <- paste0("TSLM + Fourier K=", k)
  bind_rows(
    left_join(
      score_fable_cv(
        cv_data |> model(candidate = ARIMA(precip ~ fourier(K = k) + pdq() + PDQ(0, 0, 0))),
        arima_label
      ),
      score_fable_holdout(
        train |> model(candidate = ARIMA(precip ~ fourier(K = k) + pdq() + PDQ(0, 0, 0))),
        arima_label
      ),
      by = "specification"
    ),
    left_join(
      score_fable_cv(
        cv_data |> model(candidate = TSLM(precip ~ trend() + fourier(K = k))),
        tslm_label
      ),
      score_fable_holdout(
        train |> model(candidate = TSLM(precip ~ trend() + fourier(K = k))),
        tslm_label
      ),
      by = "specification"
    )
  )
})

score_tbats_candidate <- function(use_trend, use_damped, label) {
  candidate_cv <- purrr::map_dfr(origins, function(i) {
    tr <- train |> slice(1:i)
    te <- train |> slice((i + 1):(i + h)) |> pull(precip)
    candidate <- forecast::tbats(
      ts(tr$precip, frequency = 12),
      use.box.cox = NULL,
      use.trend = use_trend,
      use.damped.trend = use_damped,
      seasonal.periods = 12
    )
    candidate_fc <- forecast::forecast(candidate, h = h)
    candidate_acc <- forecast::accuracy(candidate_fc, te, d = 0, D = 1)
    tibble(
      MASE = candidate_acc["Test set", "MASE"],
      RMSE = candidate_acc["Test set", "RMSE"],
      MAE = candidate_acc["Test set", "MAE"],
      MAPE = candidate_acc["Test set", "MAPE"]
    )
  })

  candidate <- forecast::tbats(
    train_ts,
    use.box.cox = NULL,
    use.trend = use_trend,
    use.damped.trend = use_damped,
    seasonal.periods = 12
  )
  candidate_fc <- forecast::forecast(candidate, h = h)
  candidate_acc <- forecast::accuracy(candidate_fc, test_actual, d = 0, D = 1)

  tibble(
    specification = label,
    CV_MASE = mean(candidate_cv$MASE),
    CV_RMSE = mean(candidate_cv$RMSE),
    CV_MAE = mean(candidate_cv$MAE),
    CV_MAPE = mean(candidate_cv$MAPE),
    SD_MASE = sd(candidate_cv$MASE),
    Min_MASE = min(candidate_cv$MASE),
    Max_MASE = max(candidate_cv$MASE),
    folds = nrow(candidate_cv),
    Test_MASE = candidate_acc["Test set", "MASE"],
    Test_RMSE = candidate_acc["Test set", "RMSE"],
    Test_MAE = candidate_acc["Test set", "MAE"],
    Test_MAPE = candidate_acc["Test set", "MAPE"]
  )
}

tbats_no_trend_audit <- tibble(
  specification = "TBATS no trend",
  CV_MASE = mean(cv_tbats$MASE),
  CV_RMSE = mean(cv_tbats$RMSE),
  CV_MAE = mean(cv_tbats$MAE),
  CV_MAPE = mean(cv_tbats$MAPE),
  SD_MASE = sd(cv_tbats$MASE),
  Min_MASE = min(cv_tbats$MASE),
  Max_MASE = max(cv_tbats$MASE),
  folds = nrow(cv_tbats),
  Test_MASE = acc_tbats_full["Test set", "MASE"],
  Test_RMSE = acc_tbats_full["Test set", "RMSE"],
  Test_MAE = acc_tbats_full["Test set", "MAE"],
  Test_MAPE = acc_tbats_full["Test set", "MAPE"]
)

tbats_audit <- bind_rows(
  tbats_no_trend_audit,
  score_tbats_candidate(NULL, NULL, "TBATS automatic trend")
)

model_selection_audit <- bind_rows(ets_audit, fourier_audit, tbats_audit) |>
  arrange(CV_MASE)

model_selection_family_best <- model_selection_audit |>
  mutate(family = case_when(
    grepl("^ARIMA", specification) ~ "ARIMA + Fourier",
    grepl("^TSLM", specification) ~ "TSLM + Fourier",
    grepl("^ETS", specification) ~ "ETS",
    grepl("^TBATS", specification) ~ "TBATS"
  )) |>
  group_by(family) |>
  slice_min(CV_MASE, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(CV_MASE)

cat("\n--- Candidate Specification Audit ---\n")
print(model_selection_family_best)
write.csv(
  model_selection_audit,
  "output/tables/model_selection_audit.csv",
  row.names = FALSE
)
}

if (identical(tolower(Sys.getenv("RUN_MODEL_AUDIT", "false")), "true")) {
  run_model_selection_audit()
} else if (!file.exists("output/tables/model_selection_audit.csv")) {
  message("Model audit skipped. Set RUN_MODEL_AUDIT=true to generate it.")
}
