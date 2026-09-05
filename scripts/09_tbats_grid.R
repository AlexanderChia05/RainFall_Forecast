# 09_tbats_grid.R - TBATS family grid search, screened on a
# degrees-of-freedom-corrected Ljung-Box test.
#
# Purpose: the selected TBATS configuration (no trend) fails the
# Ljung-Box residual whiteness test at lag 24 once the test is given the
# model's own degrees of freedom (forecast::modeldf() = 4, so df = 20,
# p = 0.0276). This script asks whether ANY configuration in the TBATS
# family clears p > 0.05 at lag 24 under that same corrected test, or
# whether the failure is a property of the family on this series.
#
# Two phases, matching the approach used elsewhere in this project:
#   Phase 1 - fit every configuration once on the training split and
#             screen on the corrected Ljung-Box p-values (cheap).
#   Phase 2 - run rolling-origin CV only on the survivors (expensive),
#             so overfitting is judged only for models that are
#             diagnostically admissible in the first place.
#
# Both the corrected p-value (fitdf = modeldf) and the uncorrected one
# (fitdf = 0, which is what forecast::checkresiduals() reports for a
# TBATS object) are recorded, because the pass/fail verdict depends on
# which convention is used and that fork should be visible, not hidden.
#
# Run after 01_data_pull.R. Writes output/tbats_grid_phase1.csv and,
# if any configuration survives, output/tbats_grid_phase2_cv.csv.

pkgs <- c("fpp3", "zoo", "forecast")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
library(purrr)
library(tidyr)
library(fpp3)

rain  <- readRDS("data/rain.rds")
h     <- 12
train <- rain |> filter(month <= max(month) - h)

train_ts    <- ts(train$precip, frequency = 12)
test_actual <- rain |> filter(month > max(train$month)) |> pull(precip)

acf_out_of_bounds <- function(resid, lag.max = 12) {
  r  <- na.omit(resid)
  ci <- 1.96 / sqrt(length(r))
  sum(abs(acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]) > ci)
}

# "AUTO" leaves the argument at NULL so TBATS picks by AIC.
opt <- function(x) if (identical(x, "AUTO")) NULL else as.logical(x)

# Ljung-Box that degrades gracefully: the test is undefined when the
# model uses up as many degrees of freedom as there are lags.
lb <- function(r, lag, fitdf) {
  if (fitdf >= lag) return(NA_real_)
  Box.test(r, lag = lag, type = "Ljung-Box", fitdf = fitdf)$p.value
}

# ------------------------------------------------------------ the grid
# damped trend is only meaningful when a trend state can exist, so the
# trend = FALSE branch carries a single damped setting instead of three.
grid <- bind_rows(
  expand.grid(box_cox = c("TRUE", "FALSE", "AUTO"),
              trend   = c("TRUE", "AUTO"),
              damped  = c("TRUE", "FALSE", "AUTO"),
              arma    = c("TRUE", "FALSE"),
              stringsAsFactors = FALSE),
  expand.grid(box_cox = c("TRUE", "FALSE", "AUTO"),
              trend   = "FALSE",
              damped  = "FALSE",
              arma    = c("TRUE", "FALSE"),
              stringsAsFactors = FALSE)
) |> as_tibble()

cat("Fitting", nrow(grid), "TBATS configurations...\n")

# --------------------------------------------------- phase 1: fit + LB
phase1 <- pmap_dfr(grid, function(box_cox, trend, damped, arma) {
  cat(".")
  fit <- try(forecast::tbats(train_ts,
                             use.box.cox      = opt(box_cox),
                             use.trend        = opt(trend),
                             use.damped.trend = opt(damped),
                             use.arma.errors  = as.logical(arma),
                             seasonal.periods = 12),
             silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(tibble(box_cox, trend, damped, arma, fitted = FALSE))
  }

  r   <- residuals(fit)
  dof <- forecast::modeldf(fit)
  fc  <- forecast::forecast(fit, h = h)
  acc <- forecast::accuracy(fc, test_actual)

  # Rebuild the model label the way forecast's own print method does:
  # TBATS(lambda, {p,q}, phi, {<period,k>})
  spec <- paste0(
    "TBATS(",
    if (is.null(fit$lambda)) "-" else round(fit$lambda, 3), ", {",
    length(fit$ar.coefficients), ",", length(fit$ma.coefficients), "}, ",
    if (is.null(fit$damping.parameter)) "-" else round(fit$damping.parameter, 3),
    ", {<12,", paste(fit$k.vector, collapse = ","), ">})"
  )

  tibble(
    box_cox, trend, damped, arma,
    fitted        = TRUE,
    spec          = spec,
    has_trend     = !is.null(fit$beta),
    n_ar          = length(fit$ar.coefficients),
    n_ma          = length(fit$ma.coefficients),
    k             = paste(fit$k.vector, collapse = ","),
    lambda        = if (is.null(fit$lambda)) NA_real_ else fit$lambda,
    AIC           = fit$AIC,
    model_dof     = dof,
    lb_p12        = lb(r, 12, dof),
    lb_p24        = lb(r, 24, dof),
    lb_p12_nodof  = lb(r, 12, 0),
    lb_p24_nodof  = lb(r, 24, 0),
    n_lags_out_12 = acf_out_of_bounds(r, 12),
    n_lags_out_24 = acf_out_of_bounds(r, 24),
    MASE_train    = acc["Training set", "MASE"],
    RMSE_train    = acc["Training set", "RMSE"],
    MASE_test     = acc["Test set", "MASE"],
    RMSE_test     = acc["Test set", "RMSE"]
  )
})
cat("\n")

phase1 <- phase1 |>
  mutate(passes_lb24 = !is.na(lb_p24) & lb_p24 > 0.05,
         passes_lb12 = !is.na(lb_p12) & lb_p12 > 0.05) |>
  arrange(desc(passes_lb24), desc(lb_p24))

print(phase1 |>
        select(spec, model_dof, lb_p12, lb_p24, lb_p24_nodof,
               passes_lb24, MASE_test, RMSE_test),
      n = 100)

dir.create("output", showWarnings = FALSE)
write.csv(phase1, "output/tbats_grid_phase1.csv", row.names = FALSE)

cat("\nConfigurations clearing Ljung-Box lag 24 (corrected dof):",
    sum(phase1$passes_lb24), "of", sum(phase1$fitted), "fitted\n")
cat("Configurations clearing it only without the dof correction:",
    sum(!phase1$passes_lb24 & phase1$lb_p24_nodof > 0.05, na.rm = TRUE), "\n")

# ----------------------------------------------- phase 2: CV survivors
survivors <- phase1 |> filter(passes_lb24)

if (nrow(survivors) == 0) {
  cat("\nNo configuration passes the corrected Ljung-Box test at lag 24.\n",
      "The failure is a property of the TBATS family on this series, not\n",
      "of the particular configuration that was selected.\n")
} else {
  cat("\nRunning rolling-origin CV on", nrow(survivors), "survivor(s).\n")

  origins <- seq(360, nrow(train) - h, by = 6)

  cv <- pmap_dfr(survivors |> select(box_cox, trend, damped, arma),
                 function(box_cox, trend, damped, arma) {
    folds <- map_dfr(origins, function(i) {
      tr <- train |> slice(1:i)
      te <- train |> slice((i + 1):(i + h)) |> pull(precip)
      m  <- forecast::tbats(ts(tr$precip, frequency = 12),
                            use.box.cox      = opt(box_cox),
                            use.trend        = opt(trend),
                            use.damped.trend = opt(damped),
                            use.arma.errors  = as.logical(arma),
                            seasonal.periods = 12)
      a  <- forecast::accuracy(forecast::forecast(m, h = h), te)
      tibble(MASE = a["Test set", "MASE"], RMSE = a["Test set", "RMSE"])
    })
    tibble(box_cox, trend, damped, arma,
           MASE_cv = mean(folds$MASE), RMSE_cv = mean(folds$RMSE),
           sd_MASE_cv = sd(folds$MASE), n_folds = nrow(folds))
  })

  phase2 <- survivors |>
    select(box_cox, trend, damped, arma, spec, model_dof,
           lb_p12, lb_p24, MASE_train, RMSE_train) |>
    left_join(cv, by = c("box_cox", "trend", "damped", "arma")) |>
    mutate(rmse_ratio_cv = RMSE_cv / RMSE_train,
           gap_pct_cv    = abs(MASE_cv - MASE_train) / MASE_cv,
           within_1_3x   = rmse_ratio_cv <= 1.3,
           within_10pct  = gap_pct_cv <= 0.10) |>
    arrange(rmse_ratio_cv)

  print(phase2, n = 100)
  write.csv(phase2, "output/tbats_grid_phase2_cv.csv", row.names = FALSE)
}

cat("\nDone.\n")
