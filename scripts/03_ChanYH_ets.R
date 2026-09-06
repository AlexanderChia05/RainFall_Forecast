# ETS(A,N,A) rainfall forecast

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
  snaive = SNAIVE(precip),
  ets    = ETS(precip ~ error("A") + trend("N") + season("A"))
)
fc <- fit |> forecast(h = h)

# Accuracy
accuracy_display <- bind_rows(
  fit |> accuracy() |> mutate(data_set = "Training set"),
  fc |> accuracy(rain) |> mutate(data_set = "Test set")
) |>
  transmute(
    model = recode(.model, snaive = "Seasonal naive", ets = "ETS"),
    data_set, RMSE, MAE, MAPE, MASE
  ) |>
  arrange(model, factor(data_set, levels = c("Training set", "Test set"))) |>
  mutate(across(c(RMSE, MAE, MAPE, MASE), ~ round(.x, 3)))

cat("\n--- ETS and Benchmark Accuracy ---\n")
print(accuracy_display)
fit |> select(ets) |> report()
fit |> select(ets) |> gg_tsresiduals()

# Save the fitted smoothing parameters and initial states printed by tidy().
ets_parameters <- fit |>
  select(ets) |>
  tidy() |>
  mutate(model = "ETS(A,N,A)") |>
  relocate(model)

# Residual checks
cat("\n== Ljung-Box test: ETS residuals ==\n")
print(augment(fit) |> filter(.model == "ets") |> features(.innov, ljung_box, lag = 12))
print(augment(fit) |> filter(.model == "ets") |> features(.innov, ljung_box, lag = 24))

cat("\n== Residual ACF summary ==\n")
print(augment(fit) |> filter(.model == "ets") |> as_tibble() |>
        summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
                  n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)))

# Train/test comparison
acc_train <- fit |> accuracy() |> filter(.model == "ets") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test  <- fc |> accuracy(rain) |> filter(.model == "ets") |>
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
  model(ets = ETS(precip ~ error("A") + trend("N") + season("A"))) |>
  forecast(h = h) |>
  accuracy(train, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE  = min(MASE),  max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds  = n())
print(cv_summary)

results <- tibble(
  model         = "ets_additive",
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
write.csv(ets_parameters, "output/tables/model_details/ets_parameters.csv", row.names = FALSE)

# Forecast plot
plot_from       <- yearmonth("2013 Jan")
test_actual_tbl <- rain |> as_tibble() |> filter(month > max(train$month)) |>
  transmute(month, precip)

p_fc <- fc |> filter(.model == "ets") |>
  autoplot(rain |> filter(month >= plot_from), level = c(80, 95)) +
  geom_line(data = test_actual_tbl, aes(x = month, y = precip),
            color = "red", linewidth = 0.45) +
  labs(title = "ETS(A,N,A): Forecast vs Actual", y = "mm/day", x = NULL) +
  theme_minimal()
save_fable_forecast_plot(
  fc |> filter(.model == "ets"), rain, "ETS",
  "ETS(A,N,A): Forecast vs Actual",
  "output/plots/group_summary/fc_ets.png"
)
u_ets <- augment(fit) |> filter(.model == "ets") |> as_tibble()
save_residual_diagnostic(u_ets$.innov, u_ets$month, "ETS",
                         "output/plots/group_summary/resid_ets.png")

cat("\nDone. Wrote the ETS parameter table and 2 plots.\n")
