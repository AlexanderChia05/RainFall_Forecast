source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

# EDA
rain |> autoplot(precip) + labs(title = "KL monthly precipitation (mm/day)", y = "mm/day")
rain |> gg_season(precip) + labs(title = "Seasonal plot - monsoon cycle")
rain |> ACF(precip, lag_max = 36) |> autoplot()

# Stationarity + white-noise check 
adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")

# Train/test split - last 12 months held out, no random split
h <- 12
train <- rain |> filter(month <= max(month) - h)

fit <- train |> model(
  snaive = SNAIVE(precip),
  ets    = ETS(precip ~ error("A") + trend("N") + season("A"))
)
fc <- fit |> forecast(h = h)
fc |> accuracy(rain) |> select(.model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
fc |> autoplot(rain, level = c(80, 95))

fit |> select(ets) |> report()
fit |> select(ets) |> gg_tsresiduals()

# Ljung-Box
augment(fit) |> filter(.model == "ets") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check
augment(fit) |> filter(.model == "ets") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check
acc_train <- fit |> accuracy() |> filter(.model == "ets") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "ets") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/group_summary/resid_A_ets.png",
       fit |> select(ets) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
