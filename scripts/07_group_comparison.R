# 07_group_comparison.R - group leader script. SHARED TOPIC: KL monthly
# precipitation (mm/day), NASA POWER 1981-2025 (540 obs). SDG 13.
# Combines the 4 final models (one per member, every one verified to
# pass Ljung-Box at BOTH lag=12/24 AND clear the group's RMSE-ratio
# <1.3x rule - see 03-06 headers for the family-selection evidence)
# into one train/test comparison.
#
# Member D (tbats_D) is NOT a fable model (forecast::tbats(), see
# 06_ChiaZY_tbats.R header for why it replaced STL) - it's fit and
# scored separately below, then merged into the same summary tables as
# the 3 fable members (ets_A, arima_B, tslm_C) using matching column
# names throughout (same bridging pattern as BATS in the ozone/v1
# project's group script).

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

h <- 12
train <- rain |> filter(month <= max(month) - h)

# The 3 fable-native models, one per member (see 03-05 headers for the
# diagnostics/evidence behind each pick):
#   A - ETS (trend explicitly forced off, not auto)      (03_ChanYH_ets.R)
#   B - ARIMA + Fourier(K=4), seasonal search off        (04_StephQF_arima.R)
#   C - TSLM (trend + Fourier(K=4))                      (05_HamGQ_tslm.R)
# Member D (TBATS, not fable-native) is fit separately below - see
# 06_ChiaZY_tbats.R header.
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

lb12 <- augment(fit) |> features(.innov, ljung_box, lag = 12) |>
  rename(member = .model, lb_stat_12 = lb_stat, lb_pvalue_12 = lb_pvalue)
lb24 <- augment(fit) |> features(.innov, ljung_box, lag = 24) |>
  rename(member = .model, lb_stat_24 = lb_stat, lb_pvalue_24 = lb_pvalue)

acf_check <- fit |>
  augment() |>
  as_tibble() |>
  group_by(.model) |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)) |>
  rename(member = .model)

acc_train <- fit |> accuracy() |>
  select(member = .model, MASE_train = MASE, RMSE_train = RMSE)

# --- Member D: TBATS (forecast::, not fable) -------------------------------
# Single holdout fit + score, matching the fable members' train/test split
# exactly (same h=12, same train cutoff). use.trend=FALSE per
# 06_ChiaZY_tbats.R's evidence (tbats_notrend beat tbats_auto on CV).
train_ts   <- ts(train$precip, frequency = 12)
fit_tbats  <- forecast::tbats(train_ts, use.box.cox = NULL, use.trend = FALSE,
                               use.damped.trend = FALSE, seasonal.periods = 12)
fc_tbats   <- forecast::forecast(fit_tbats, h = h)
test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)
acc_tbats_full <- forecast::accuracy(fc_tbats, test_actual)
resid_tbats <- residuals(fit_tbats)
lb_tbats12  <- Box.test(resid_tbats, lag = 12, type = "Ljung-Box")
lb_tbats24  <- Box.test(resid_tbats, lag = 24, type = "Ljung-Box")

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

# Rolling-origin CV (.init=360, .step=6). stretch_tsibble() generates
# windows up to n=540 inclusive (31 total, k=0..30, origin=360+6k) - but
# only k=0..28 (origin<=528) have a FULL 12-month actual future to score
# against. k=29 (origin=534) only has 6 of 12 real months (535-540) and
# still produces a non-NA MASE/RMSE from that partial horizon - filtering
# on !is.na(MASE) alone lets this partial fold sneak into the average.
# k=30 (origin=540) has 0 real future months and correctly comes back NA.
# Explicitly cap at n_folds_clean=29 (.id <= 29, i.e. origin<=528) so
# every fold here has a genuine complete 12-month horizon, matching
# TBATS's manual loop exactly (both now compare the identical 29 origins,
# 360 through 528) - not an approximate match, an exact one.
n_folds_clean <- length(seq(360, nrow(rain) - h, by = 6))  # 29
origins <- seq(360, nrow(rain) - h, by = 6)

set.seed(2026)
cv_fits <- rain |>
  stretch_tsibble(.init = 360, .step = 6) |>
  model(
    snaive  = SNAIVE(precip),
    ets_A   = ETS(precip ~ error("A") + trend("N") + season("A")),
    arima_B = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
    tslm_C  = TSLM(precip ~ trend() + fourier(K = 4))
  )

cv_acc <- cv_fits |> forecast(h = 12) |> accuracy(rain, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  filter(!is.na(MASE), .id <= n_folds_clean) |>
  group_by(member = .model) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE = min(MASE), max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds = n())

# TBATS CV: refit at each of the same origins (slower than the fable
# picks - each fold is a fresh forecast::tbats() call).
cv_tbats <- map_dfr(origins, function(i) {
  tr  <- rain |> slice(1:i)
  te  <- rain |> slice((i + 1):(i + h)) |> pull(precip)
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

# Overfitting check (MUST) - both the group's original MASE-gap
# convention and the RMSE-ratio<1.3x rule, both computed against the CV
# mean (authoritative), not the single holdout (reference only, kept
# alongside for transparency).
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
