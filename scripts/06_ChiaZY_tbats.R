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

train_ts <- ts(train$precip, frequency = 12)

fit_tbats <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                              use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats  <- forecast::forecast(fit_tbats, h = h)
print(fit_tbats)

test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats <- forecast::accuracy(fc_tbats, test_actual)
print(acc_tbats)

resid_tbats <- residuals(fit_tbats)

# Ljung-Box
Box.test(resid_tbats, lag = 12, type = "Ljung-Box")
Box.test(resid_tbats, lag = 24, type = "Ljung-Box")

# ACF-in-bounds check
acf_out_of_bounds(resid_tbats, lag.max = 12)
acf_out_of_bounds(resid_tbats, lag.max = 24)

# Overfitting check
mase_train <- acc_tbats["Training set", "MASE"]
mase_test  <- acc_tbats["Test set", "MASE"]
rmse_train <- acc_tbats["Training set", "RMSE"]
rmse_test  <- acc_tbats["Test set", "RMSE"]
cat("rmse_ratio_holdout:", rmse_test / rmse_train,
    " mase_gap_pct_holdout:", abs(mase_test - mase_train) / mase_test, "\n")

dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
png("output/plots/group_summary/resid_D_tbats.png", width = 800, height = 600, res = 150)
forecast::checkresiduals(fit_tbats)
dev.off()
