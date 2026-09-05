# TBATS Codebase Review

## Review conclusion

The handoff document is only partly correct and should not be applied as a complete set of fixes.

The MASE scaling problem is real in the `main` branch. However, it has already been fixed in the current `steph` branch. The proposed Ljung-Box degrees-of-freedom correction is more questionable because it differs from the default approach used by `forecast::checkresiduals()` for TBATS. As a result, the recommendation to disable the Box-Cox transformation should not be accepted without first choosing and justifying the Ljung-Box convention.

Under the convention currently used by the `forecast` package, the latest TBATS model passes the reported residual tests at the 5% level and achieves the lowest cross-validation and test errors. The current conclusion that TBATS is the best overall forecasting model is therefore supported by the latest output.

## Implementation update (2026-09-06)

The reporting and reproducibility changes from this review have now been implemented and executed:

- All accuracy calculations retain seasonal MASE scaling with `d = 0` and `D = 1`.
- TBATS remains automatic Box-Cox with the selected lambda of approximately 0.596.
- `output/tables/model_specifications.csv` records the comparable specification and Ljung-Box convention for all models.
- Model-specific parameter and specification CSV files are produced for ETS, ARIMA, TSLM and TBATS.
- `output/tables/data_summary.csv` records the data source, period, location and missing-value counts.
- `output/tables/package_versions.csv` records the R and forecasting-package versions used for the final run.
- The EDA script saves the seasonal, subseries, ACF and PACF figures instead of leaving them only in an interactive plotting device.
- The full comparison and all four individual model scripts completed successfully. The latest ranking and diagnostic conclusions are unchanged.

## 1 MASE scaling

### Finding

The handoff correctly identifies a MASE comparability problem in the earlier `main` branch.

When test observations are supplied to `forecast::accuracy()`, the function uses the frequency of the supplied `x` object to choose its default MASE denominator. In this codebase, `test_actual` and `te` are created with `pull()`. They are ordinary numeric vectors and have a frequency of one. If `d` and `D` are omitted, `forecast::accuracy()` treats them as non-seasonal data and uses a lag-one naive denominator.

For this monthly series, the intended denominator is based on a seasonal difference at lag 12. The correct explicit settings are:

```r
d = 0
D = 1
```

Official source:

- https://github.com/robjhyndman/forecast/blob/master/R/errors.R

### Current status

The current `steph` branch already supplies these arguments in all four required locations:

- `scripts/06_ChiaZY_tbats.R`, lines 108 to 113
- `scripts/06_ChiaZY_tbats.R`, lines 162 to 167
- `scripts/07_group_comparison.R`, lines 53 to 58
- `scripts/07_group_comparison.R`, lines 123 to 128

The current output also confirms that seasonal scaling is being used:

```text
TBATS test MAE / TBATS test MASE
= 1.983726 / 0.770576
= approximately 2.574
```

The resulting denominator matches the seasonal naive training MAE of approximately 2.574.

### Verdict

The MASE bug exists in the earlier `main` version but is already fixed in the current `steph` code. No additional change is required for Issue 3.

## 2 Ljung-Box degrees of freedom

### Handoff claim

The handoff proposes the following correction:

```r
tbats_dof <- forecast::modeldf(fit_tbats)

Box.test(
  resid_tbats,
  lag = 24,
  type = "Ljung-Box",
  fitdf = tbats_dof
)
```

For the automatic Box-Cox model, `modeldf(fit_tbats)` returns four. This changes the lag 24 p-value from approximately .0886 to .0276.

### Review finding

This should not automatically be treated as the correct TBATS diagnostic.

The official implementation of `forecast::checkresiduals()` applies a model degrees-of-freedom adjustment automatically for ARIMA models. For other models, including TBATS, it uses `fitdf = 0` for the Ljung-Box test.

Official source:

- https://github.com/robjhyndman/forecast/blob/master/R/checkresiduals.R
- https://pkg.robjhyndman.com/forecast/reference/checkresiduals.html

The official R documentation for `Box.test()` specifically describes the usual adjustment for residuals from an ARMA(p,q) model as:

```text
fitdf = p + q
```

Official documentation:

- https://stat.ethz.ch/R-manual/R-devel/library/stats/help/Box.test.html

The fitted TBATS model is reported as:

```text
TBATS(0.596, {0,0}, -, {<12,3>})
```

Its ARMA order is `{0,0}`. An ARMA-based degrees-of-freedom correction would therefore give `fitdf = 0`.

Using `forecast::modeldf(fit_tbats) = 4` subtracts Box-Cox and state-space parameters in addition to any ARMA terms. This is more conservative than both `forecast::checkresiduals()` and the ARMA guidance in the base R documentation. It may be reported as a sensitivity check, but it is not clearly the only correct test.

### Verdict

The absence of `fitdf = forecast::modeldf(fit_tbats)` in the current code is not a confirmed bug.

For this assignment, the recommended primary convention is the package-standard TBATS diagnostic with `fitdf = 0`. Under this convention, the lag 24 p-value is approximately .089. The correct interpretation is that the test does not reject the null hypothesis of no residual autocorrelation at the 5% level. It does not prove that the residuals are fully independent.

## 3 Box-Cox model selection

### Current specification

Both current scripts fit TBATS using:

```r
forecast::tbats(
  train_ts,
  use.box.cox = NULL,
  use.trend = FALSE,
  use.damped.trend = FALSE,
  seasonal.periods = 12
)
```

`use.box.cox = NULL` allows the algorithm to compare transformed and untransformed fits using AIC. The current fitted model selects a Box-Cox parameter of approximately 0.596.

### Current results

| Measure | Result |
|---|---:|
| Training RMSE | 2.449 |
| Training MASE | 0.730 |
| CV RMSE | 2.495 |
| CV MASE | 0.783 |
| Test RMSE | 2.335 |
| Test MAE | 1.984 |
| Test MAPE | 35.702 |
| Test MASE | 0.771 |
| Ljung-Box p-value at lag 12 | .524 |
| Ljung-Box p-value at lag 24 | .089 |

TBATS has the lowest CV RMSE and CV MASE among the four main models. It also has the lowest RMSE, MAE, MAPE and MASE on the 2025 test set.

### Assessment of the handoff recommendation

The `main` branch changed `use.box.cox` to `FALSE` because the automatic Box-Cox model failed when `modeldf(fit_tbats) = 4` was subtracted from the Ljung-Box degrees of freedom. Since that degrees-of-freedom convention is disputed, the argument for forcing `use.box.cox = FALSE` is also not decisive.

Under the package-standard `fitdf = 0` diagnostic, the automatic Box-Cox model has p = .089 and passes at the 5% level. It also provides better forecast accuracy. The current `use.box.cox = NULL` setting is therefore reasonable.

### Verdict

Retain `use.box.cox = NULL` unless the lecturer specifically requires all TBATS state parameters to be deducted from the Ljung-Box degrees of freedom.

Do not combine the following results because they refer to different model and evaluation conventions:

- Old no-Box-Cox CV MASE of 0.701
- Old corrected Ljung-Box p-value of .0587
- Current automatic Box-Cox CV MASE of 0.783
- Current uncorrected Ljung-Box p-value of .089

## 4 Branch relationship

The handoff describes `steph` as a parallel rewrite that branched before the main fixes. Git history does not support that description.

The current history is:

```text
09d8e61  steph  enhance for graphs and tbats command
8136d8e  main   match group comparison TBATS specification
cc3a1fc         switch TBATS to Box-Cox-off
```

The merge base of `main` and `steph` is commit `8136d8e`. Therefore, `steph` is a direct one-commit descendant of `main`. Commit `09d8e61` deliberately changed the TBATS configuration back to automatic Box-Cox, added the explicit seasonal MASE arguments and changed the output naming scheme.

### Verdict

The recommendation to discard `steph` as an unrelated rewrite is not justified. The individual changes in commit `09d8e61` should be evaluated separately.

## 5 Output consistency

The current individual TBATS script and group comparison script use the same settings:

- `use.box.cox = NULL`
- `use.trend = FALSE`
- `use.damped.trend = FALSE`
- `seasonal.periods = 12`
- MASE with `d = 0` and `D = 1`
- Ljung-Box with `fitdf = 0`

The following latest files are internally consistent with that specification:

- `output/tables/model_results.csv`
- `output/tables/model_specifications.csv`
- `output/tables/model_details/tbats_parameters.csv`
- `output/plots/group_summary/fc_tbats.png`
- `output/plots/group_summary/resid_tbats.png`

The current model ranking is therefore usable as long as the package-standard Ljung-Box convention is retained.

## 6 Additional omissions

### Model specification outputs — resolved

The scripts now save complete model specifications rather than relying only on console output. TBATS records its selected Box-Cox lambda, ARMA order, harmonics, trend settings, MASE convention and Ljung-Box convention. ARIMA records its selected orders and Fourier setting. ETS and TSLM parameter tables are also saved. These files can now be used directly when preparing the individual-report appendices.

### Bias adjustment under Box-Cox

When a Box-Cox transformation is used, `forecast::forecast()` distinguishes between a back-transformed median forecast and a bias-adjusted mean forecast. The default behaviour normally produces a median when bias adjustment is disabled.

Official documentation:

- https://search.r-project.org/CRAN/refmans/forecast/html/forecast.bats.html

This is not a critical error for the current assignment. However, RMSE is theoretically associated with mean forecasts. A sensitivity check could compare the current result with:

```r
forecast::forecast(fit_tbats, h = h, biasadj = TRUE)
```

If this option is tested, it must be applied consistently to the full fit and every CV fold. All TBATS accuracy outputs and forecast figures must then be regenerated.

### Reproducibility

The standalone member scripts download NASA POWER data again, while the group comparison uses `data/rain.rds`. NASA may revise historical grid values. Running scripts on different dates could therefore create small differences between member and group outputs.

For the final submission, all scripts should use the same saved dataset or the dataset should be downloaded once immediately before all scripts are executed.

Package installation also occurs automatically inside the scripts. Package updates could change model selection or output. The final comparison script now records the versions used in `output/tables/package_versions.csv`.

## Recommended final decision

1. Keep `d = 0` and `D = 1` in all TBATS accuracy calculations.
2. Keep `use.box.cox = NULL` in both the individual and group scripts.
3. Use `fitdf = 0` as the primary TBATS Ljung-Box convention because it matches `forecast::checkresiduals()` and the fitted model has ARMA order `{0,0}`.
4. Interpret p = .089 as a diagnostic pass at the 5% level, not proof of independent residuals.
5. Retain the current conclusion that TBATS is the best overall model based on CV and 2025 test accuracy.
6. Do not use the old grid-search values or the old no-Box-Cox results in the reports.
7. Use the newly generated model specification and parameter CSV files in the individual-report appendices.
8. Consider bias adjustment only as a sensitivity check. It is not required for the present beginner-level report.

## Files reviewed

- `D:/RDS/Y2S1/BMMS2094 Statistics for Data Science/asm/CODEX_HANDOFF_TBATS_ISSUES.md`
- `scripts/06_ChiaZY_tbats.R`
- `scripts/07_group_comparison.R`
- `scripts/00_setup.R`
- `scripts/01_data_pull.R`
- Current CSV outputs under `output/`
- Current forecast and residual plots under `output/plots/group_summary/`
