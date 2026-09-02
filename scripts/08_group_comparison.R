# 08_group_comparison.R - group leader script. SHARED TOPIC: KL monthly
# precipitation (mm/day), NASA POWER 1981-2025 (540 obs). SDG 13.
# Combines the 4 final models (one per member, every one verified via
# 03_family_grid_search.R to pass Ljung-Box at BOTH lag=12/24 AND clear
# the group's RMSE-ratio<1.3x rule - see 02-05 headers for the
# family-selection evidence) into one train/test comparison.
#
# All 4 families are fable-native this project (no forecast::-package
# bridging needed, unlike BATS in the ozone/v1 project) - the whole
# comparison fits in one mable, no special-casing required.

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

h <- 12
train <- rain |> filter(month <= max(month) - h)

# The 4 final models, one per member (see 02-05 headers for the
# diagnostics/evidence behind each pick):
#   A - ETS (auto-selected, converges to no-trend)      (04_ChanYH_ets.R)
#   B - ARIMA + Fourier(K=4), seasonal search off        (05_StephQF_arima.R)
#   C - TSLM (trend + seasonal dummies)                  (06_HamGQ_tslm.R)
#   D - STL(robust) + ARIMA(remainder)                   (07_ChiaZY_stl.R)
set.seed(2026)
fit <- train |> model(
  snaive      = SNAIVE(precip),
  ets_A       = ETS(precip),
  arima_B     = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
  tslm_C      = TSLM(precip ~ trend() + season()),
  stl_D       = decomposition_model(STL(precip, robust = TRUE),
                                     ARIMA(season_adjust ~ pdq()))
)

fc <- fit |> forecast(h = h)
acc_test <- fc |> accuracy(rain) |>
  select(member = .model, MASE, RMSE, MAE, MAPE) |> arrange(MASE)
print(acc_test)

# Ljung-Box p-value per model - residuals must look random (p > 0.05).
lb12 <- augment(fit) |> features(.innov, ljung_box, lag = 12) |>
  rename(member = .model, lb_stat_12 = lb_stat, lb_pvalue_12 = lb_pvalue)
lb24 <- augment(fit) |> features(.innov, ljung_box, lag = 24) |>
  rename(member = .model, lb_stat_24 = lb_stat, lb_pvalue_24 = lb_pvalue)
print(lb12)
print(lb24)

# ACF-in-bounds check (MUST)
acf_check <- fit |>
  augment() |>
  as_tibble() |>
  group_by(.model) |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24)) |>
  rename(member = .model)
print(acf_check)

acc_train <- fit |> accuracy() |>
  select(member = .model, MASE_train = MASE, RMSE_train = RMSE)

# Rolling-origin CV (.init=360, .step=6 -> ~29 folds). Denser than the
# ozone project's 6-fold convention from the start (not a retrofit) -
# all 4 families here are comparatively fast to fit (no NNAR/BATS), so
# the extra folds cost little and give a much more stable ratio/gap
# estimate (same lesson the ozone project learned the hard way in
# 11_denser_cv.R).
set.seed(2026)
cv_fits <- rain |>
  stretch_tsibble(.init = 360, .step = 6) |>
  model(
    snaive  = SNAIVE(precip),
    ets_A   = ETS(precip),
    arima_B = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
    tslm_C  = TSLM(precip ~ trend() + season()),
    stl_D   = decomposition_model(STL(precip, robust = TRUE),
                                   ARIMA(season_adjust ~ pdq()))
  )

cv_acc <- cv_fits |> forecast(h = 12) |> accuracy(rain, by = c(".model", ".id"))

cv_summary <- cv_acc |>
  filter(!is.na(MASE)) |>
  group_by(member = .model) |>
  summarise(mean_MASE = mean(MASE), sd_MASE = sd(MASE),
            min_MASE = min(MASE), max_MASE = max(MASE),
            mean_RMSE = mean(RMSE), n_folds = n()) |>
  arrange(mean_MASE)
print(cv_summary)
write.csv(cv_summary, "output/model_comparison_cv_summary.csv", row.names = FALSE)

# Overfitting check (MUST) - both the group's original MASE-gap
# convention (carried over from the ozone project) and the newer
# RMSE-ratio<1.3x rule, both computed against the CV mean (authoritative)
# not the single holdout (reference only, kept alongside for
# transparency).
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
