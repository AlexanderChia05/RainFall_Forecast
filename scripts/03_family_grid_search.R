# 03_family_grid_search.R - group leader script. Systematic grid search
# across ETS, ARIMA, STL, and TSLM families for KL monthly precipitation
# (SDG 13). Same two-phase design proven across v1 (ozone, 09_...) and
# v2 (trade, 02_...):
#   Phase 1 (fast, single holdout): fit every variant once, get Ljung-Box
#     p-value at BOTH lag=12 and lag=24.
#   Phase 2 (rolling CV): ONLY for variants passing BOTH lags in phase 1
#     - RMSE_test_cv/RMSE_train ratio (1.3x rule) + MASE gap percentage.
#
# Differences from the trade grid search (lessons applied, not repeated
# blind):
#   - No regime dummies needed here (no GST/SST-equivalent structural
#     break for rainfall) - simpler formulas, no forecast(new_data=...)
#     workaround required, plain forecast(h=h) works for every variant.
#   - fourier() variants all carry PDQ(0,0,0) explicitly from the start
#     (the trade run's arima_fourierK2/K4 bug - fourier vs ARIMA's own
#     seasonal search competing for the same signal, exact-singular
#     design matrix) - fixed proactively instead of hit-and-debug.
#   - trend_strength was measured weak (0.181, see 02_eda_stationarity.R)
#     - added trend("N")/season()-only variants per family rather than
#     assuming a trend term is needed everywhere.
#   - No NNAR in this family set (not one of the 4 chosen families) and
#     no BATS either - all 4 families are fable-native and comparatively
#     fast to fit, so phase 2 CV uses MORE folds (.step=6, .init=360 ->
#     ~29 folds) since compute budget allows it.

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")
h <- 12
train <- rain |> filter(month <= max(month) - h)

# ============================================================
# PHASE 1: fit every variant once (single holdout)
# ============================================================
family_specs <- list(
  # --- ETS family ---
  ets_auto        = ETS(precip),
  ets_damped      = ETS(precip ~ error("A") + trend("Ad") + season("A")),
  ets_multseason  = ETS(precip ~ error("M") + trend("N") + season("M")),
  ets_notrend     = ETS(precip ~ error("A") + trend("N") + season("A")),
  # --- ARIMA family ---
  arima_auto         = ARIMA(precip),
  arima_fourierK2    = ARIMA(precip ~ fourier(K = 2) + pdq() + PDQ(0, 0, 0)),
  arima_fourierK4    = ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0, 0, 0)),
  arima_fourierK6    = ARIMA(precip ~ fourier(K = 6) + pdq() + PDQ(0, 0, 0)),
  arima_seasonalonly = ARIMA(precip ~ pdq(0, 0, 0) + PDQ()),
  # --- STL family ---
  stl_arima     = decomposition_model(STL(precip, robust = TRUE), ARIMA(season_adjust ~ pdq())),
  stl_ets       = decomposition_model(STL(precip, robust = TRUE), ETS(season_adjust)),
  stl_notrobust = decomposition_model(STL(precip, robust = FALSE), ARIMA(season_adjust ~ pdq())),
  stl_rw        = decomposition_model(STL(precip, robust = TRUE), RW(season_adjust)),
  # --- TSLM family ---
  tslm_trendseason  = TSLM(precip ~ trend() + season()),
  tslm_seasononly   = TSLM(precip ~ season()),
  tslm_fourier      = TSLM(precip ~ trend() + fourier(K = 4)),
  tslm_fourier_nt   = TSLM(precip ~ fourier(K = 6))
)

fit_grid <- train |> model(!!!family_specs)

fc_grid   <- fit_grid |> forecast(h = h)
acc_grid  <- fc_grid |> accuracy(rain) |> select(variant = .model, MASE, RMSE, MAE, MAPE)
train_acc <- fit_grid |> accuracy() |> select(variant = .model, MASE_train = MASE, RMSE_train = RMSE)
lb12_grid <- augment(fit_grid) |> features(.innov, ljung_box, lag = 12) |>
  rename(variant = .model, lb_stat_12 = lb_stat, lb_pvalue_12 = lb_pvalue)
lb24_grid <- augment(fit_grid) |> features(.innov, ljung_box, lag = 24) |>
  rename(variant = .model, lb_stat_24 = lb_stat, lb_pvalue_24 = lb_pvalue)

ets_names  <- c("ets_auto","ets_damped","ets_multseason","ets_notrend")
arima_names<- c("arima_auto","arima_fourierK2","arima_fourierK4","arima_fourierK6","arima_seasonalonly")
stl_names  <- c("stl_arima","stl_ets","stl_notrobust","stl_rw")
tslm_names <- c("tslm_trendseason","tslm_seasononly","tslm_fourier","tslm_fourier_nt")

phase1_all <- acc_grid |> left_join(train_acc, by = "variant") |>
  left_join(lb12_grid |> select(variant, lb_stat_12, lb_pvalue_12), by = "variant") |>
  left_join(lb24_grid |> select(variant, lb_stat_24, lb_pvalue_24), by = "variant") |>
  mutate(family = case_when(variant %in% ets_names   ~ "ETS",
                             variant %in% arima_names ~ "ARIMA",
                             variant %in% stl_names   ~ "STL",
                             variant %in% tslm_names  ~ "TSLM",
                             TRUE ~ "other"),
         passes_both_lags = lb_pvalue_12 > 0.05 & lb_pvalue_24 > 0.05) |>
  select(family, variant, MASE, RMSE, MAE, MAPE, MASE_train, RMSE_train,
         lb_stat_12, lb_pvalue_12, lb_stat_24, lb_pvalue_24, passes_both_lags) |>
  arrange(family, desc(passes_both_lags), MASE)

print(phase1_all, n = Inf)
write.csv(phase1_all, "output/family_grid_phase1_lb_results.csv", row.names = FALSE)

passers <- phase1_all |> filter(passes_both_lags) |> pull(variant)
cat("\n== PHASE 1 DONE. Variants passing BOTH LB lags (going to phase 2 CV):",
    if (length(passers)) paste(passers, collapse = ", ") else "NONE", "==\n")

# ============================================================
# PHASE 2: rolling CV (.init=360, .step=6, ~29 folds), phase-1 passers only
# ============================================================
if (length(passers) > 0) {
  set.seed(2026)
  cv_fits <- rain |> stretch_tsibble(.init = 360, .step = 6) |>
    model(!!!family_specs[passers])
  cv_acc <- cv_fits |> forecast(h = 12) |> accuracy(rain, by = c(".model", ".id"))
  cv_summary <- cv_acc |> filter(!is.na(MASE)) |>
    group_by(variant = .model) |>
    summarise(mean_MASE_cv = mean(MASE), sd_MASE_cv = sd(MASE),
              mean_RMSE_cv = mean(RMSE), n_folds = n())

  final_verdict <- phase1_all |> filter(variant %in% passers) |>
    left_join(cv_summary, by = "variant") |>
    mutate(rmse_ratio    = mean_RMSE_cv / RMSE_train,
           mase_gap_pct  = abs(mean_MASE_cv - MASE_train) / mean_MASE_cv,
           fit_flag      = case_when(rmse_ratio > 1.3     ~ "Overfitting",
                                      rmse_ratio < 1 / 1.3 ~ "Underfitting",
                                      TRUE                 ~ "Acceptable")) |>
    select(family, variant, MASE_train, RMSE_train, mean_MASE_cv, mean_RMSE_cv, n_folds,
           mase_gap_pct, rmse_ratio, fit_flag, lb_pvalue_12, lb_pvalue_24) |>
    arrange(rmse_ratio)

  print(final_verdict, n = Inf)
  write.csv(final_verdict, "output/family_grid_phase2_rmse_ratio.csv", row.names = FALSE)

  cat("\n== Variants with RMSE ratio WITHIN 1.3x (both LB lags also pass): ==\n")
  print(final_verdict |> filter(fit_flag == "Acceptable"))
} else {
  cat("No variant passed both LB lags in phase 1 - phase 2 skipped, nothing to CV.\n")
}
