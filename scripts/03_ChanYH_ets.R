# 03_ChanYH_ets.R - SHARED TOPIC: KL monthly mean precipitation rate
# (mm/day), NASA POWER, 1981-2025 (540 obs). SDG 13 (Climate Action)
# primary, SDG 11/6 extensions - see 01_data_pull.R header. Non-white-
# noise (Ljung-Box on raw series p < 2.2e-16, both lags), weak trend
# (trend_strength=0.181), moderate seasonal strength (0.469, peak
# November - NE monsoon). Member A model family: ETS (exponential
# smoothing / Holt-Winters).
#
# Model pick: ETS(error("A") + trend("N") + season("A")) - trend
# EXPLICITLY forced off, not left to auto-selection. Verified via a
# 17-variant grid across ETS/ARIMA/STL/TSLM (grid search script since
# removed once the search phase was done, numbers held here for the
# record): ets_auto (unconstrained) and this ets_notrend spec landed on
# the SAME fit byte-for-byte on this data snapshot - MASE_train=0.738,
# RMSE_train=2.44, mean_MASE_cv=0.805, mean_RMSE_cv=2.54, ratio=1.04,
# gap=8.3%, p12=0.522, p24=0.122 for both (final numbers per the 29-
# fold-clean CV, see 07_group_comparison.R - a slightly earlier, pre-
# fold-parity-fix run of this same spec had shown 0.813/2.55/1.05/9.2%;
# the fold-count fix nudged these numbers slightly, not the model
# itself). Forced explicitly (rather
# than relying on ETS() to keep re-deriving the same answer) so the
# pick doesn't silently drift to a different structure if the series is
# re-pulled with a few more months and auto-selection's AICc comparison
# tips the other way - the weak trend_strength (0.181, measured in
# 02_eda_stationarity.R) is why "no trend" was right in the first
# place, not a manual simplification imposed after the fact.
#
# Why this family works here (unlike the ozone/trade topics): ETS has
# no ARMA-error term or external-regressor slot, which is exactly why
# it failed on both prior topics (ozone's QBO-driven residual
# autocorrelation, trade's business-cycle momentum - neither is a clean
# deterministic trend+season shape ETS can represent). Rainfall's
# seasonality is a genuinely mechanical, near-fixed monsoon cycle with
# no such confound, so ETS's trend+season decomposition is sufficient
# on its own here.
#
# Diagnostics (single holdout, from the grid search): p=0.522 (lag12),
# p=0.122 (lag24) - clears both conventions with real margin, unlike
# ozone's NNAR (thin 0.0594) or any ETS variant tried on ozone/trade
# (all rejected, best case ~0.00087). 29-fold rolling CV (.init=360,
# .step=6, capped to exclude the one partial-horizon fold - see
# 07_group_comparison.R): RMSE ratio = 1.04 (within the 1.3x rule),
# MASE gap = 8.3% (within the group's own 10% target - a first across
# all 3 topics attempted this project).

source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")

# EDA
rain |> autoplot(precip) + labs(title = "KL monthly precipitation (mm/day)", y = "mm/day")
rain |> gg_season(precip) + labs(title = "Seasonal plot - monsoon cycle")
rain |> ACF(precip, lag_max = 36) |> autoplot()

# Stationarity + white-noise check (see 02_eda_stationarity.R for the
# full ADF/KPSS discussion - they disagree, explained there, not
# repeated here)
adf.test(rain$precip)
kpss.test(rain$precip)
Box.test(rain$precip, lag = 12, type = "Ljung-Box")  # want p < 0.05 -> not white noise

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

# Ljung-Box - residuals must look random (p > 0.05), both conventions
augment(fit) |> filter(.model == "ets") |> features(.innov, ljung_box, lag = 12)
augment(fit) |> filter(.model == "ets") |> features(.innov, ljung_box, lag = 24)

# ACF-in-bounds check (MUST) - count residual ACF lags outside +-1.96/sqrt(n)
augment(fit) |> filter(.model == "ets") |> as_tibble() |>
  summarise(n_lags_out_12 = acf_out_of_bounds(.innov, lag.max = 12),
            n_lags_out_24 = acf_out_of_bounds(.innov, lag.max = 24))

# Overfitting check (single holdout, reference only - see
# 07_group_comparison.R for the authoritative CV-based ratio/gap)
acc_train <- fit |> accuracy() |> filter(.model == "ets") |>
  select(.model, MASE_train = MASE, RMSE_train = RMSE)
acc_test <- fc |> accuracy(rain) |> filter(.model == "ets") |>
  select(.model, MASE_test = MASE, RMSE_test = RMSE)
acc_train |> left_join(acc_test, by = ".model") |>
  mutate(rmse_ratio = RMSE_test / RMSE_train,
         mase_gap_pct = abs(MASE_test - MASE_train) / MASE_test)

dir.create("output/plots/ChanYH_ets", recursive = TRUE, showWarnings = FALSE)
ggsave("output/plots/ChanYH_ets/resid_ets.png",
       fit |> select(ets) |> gg_tsresiduals(), width = 8, height = 6, dpi = 150)
