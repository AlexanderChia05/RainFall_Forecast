source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

rain |> autoplot(precip)
rain |> gg_season(precip)
rain |> ACF(precip, lag_max = 36) |> autoplot()

adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")

h <- 12
train <- rain |> filter(month <= max(month) - h)

fit <- train |> model(
  snaive = SNAIVE(precip),
  tslm   = TSLM(precip ~ trend() + fourier(K = 4))
)
fc <- fit |> forecast(h = h)
fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
fc |> autoplot(rain, level = c(80, 95))

fit |> select(tslm) |> report()
fit |> select(tslm) |> gg_tsresiduals()

# Ljung-Box
augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 12)
augment(fit) |> filter(.model == "tslm") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check
augment(fit) |> filter(.model == "tslm") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check
acc_train <- fit |> accuracy() |> filter(.model == "tslm") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "tslm") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/group_summary/resid_C_tslm.png",
       fit |> select(tslm) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
