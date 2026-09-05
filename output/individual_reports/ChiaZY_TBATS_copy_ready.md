# Forecasting Monthly Mean Precipitation in Kuala Lumpur Using TBATS

## I METHODOLOGY

### A Data preparation and model motivation

The analysis used monthly mean precipitation rates for Kuala Lumpur at 3.1390 N and 101.6869 E. The data were downloaded from NASA POWER and were measured in millimetres per day (NASA Langley Research Center, n.d.). The observation period ran from January 1981 to December 2025. It provided 540 monthly values. Values coded as -999 were converted to missing values. Linear interpolation was included in the preprocessing procedure. No missing values remained in the final series.

The data were divided in time order. The 528 observations ending in December 2024 formed the training set. The 12 observations in 2025 formed an independent test set. Keeping the last year outside model fitting gave a direct test of how the model performs on a future period.

STL decomposition produced a seasonal strength of 0.469 and a trend strength of 0.181. Figure A1 shows regular annual movement together with changes in seasonal amplitude across the record. The raw series also failed the Ljung-Box white-noise test at lag 24 with a statistic of 233.12 and p < .001. A method that can represent annual seasonality and allow the variance structure to adjust was therefore reasonable.

### B TBATS specification and justification

TBATS represents trigonometric seasonality and may include a Box-Cox transformation, ARMA errors and trend components. It was developed for seasonal series that require a flexible state-space form (De Livera et al., 2011). The model in this study was fitted with a seasonal period of 12. Trend and damped trend were switched off because the measured trend strength was weak. The Box-Cox choice was left to the estimation procedure. ARMA error terms and the number of trigonometric components were also determined by the fitting algorithm. The implementation used the forecast package in R (Hyndman, 2024).

### C Evaluation procedure

Forecast accuracy was evaluated with RMSE, MAE, MAPE and MASE. MASE was treated as the main scale-free measure. A value below one means that the error is smaller than the in-sample seasonal naive scaling error (Hyndman & Athanasopoulos, 2021). RMSE gives more weight to large errors. MAE gives the average absolute error in mm/day. MAPE was included but was not used alone because percentage errors can be affected by months with low rainfall.

The seasonal naive method was used as a benchmark. Rolling-origin cross-validation began with 360 months. The training window increased by six months after each origin. Every fit forecasted the following 12 months. This generated 27 complete folds and stopped before the 2025 test year. For TBATS accuracy, the seasonal MASE scale was specified with d = 0 and D = 1.

## II DATA ANALYSIS RESULTS AND DISCUSSION

TBATS obtained a training RMSE of 2.449 and a training MASE of 0.730. Its mean CV RMSE was 2.495 and its mean CV MASE was 0.783. The CV RMSE was only 1.019 times the training RMSE. The MASE difference was 6.84 percent when measured against the CV value. This was the smallest training-to-CV difference among the four main models. The result suggests that TBATS did not gain its accuracy by fitting the training sample too closely.

The seasonal naive benchmark gave a CV RMSE of 3.179 and a CV MASE of 0.985. TBATS reduced CV RMSE by about 21.5 percent and CV MASE by about 20.5 percent. It achieved the lowest CV RMSE and CV MASE among TBATS, ARIMA with Fourier terms, ETS and TSLM with Fourier terms.

The independent 2025 results supported the same conclusion. TBATS recorded an RMSE of 2.335, an MAE of 1.984, a MAPE of 35.702 percent and a MASE of 0.771. All four values were the lowest in the final model comparison. The seasonal naive test RMSE was 2.694 and its test MASE was 0.942. TBATS reduced test RMSE by 13.3 percent and test MASE by 18.2 percent. Test MAE fell by 18.2 percent. Test MAPE fell by 17.3 percent. This agreement between CV and the final test year provides stronger evidence than training accuracy alone.

Figure A2 shows that the forecast follows the broad annual direction but remains smoother than the actual monthly rainfall. The actual series contains sharp changes that no model reproduced fully. All 12 observed test values appear within the plotted 95 percent interval. The intervals therefore covered the selected test year while still reflecting substantial rainfall uncertainty.

The Ljung-Box p-values were .524 at lag 12 and .089 at lag 24. Under the .05 rule used in the project, the null hypothesis of no residual autocorrelation was not rejected. The model therefore passed both reported Ljung-Box checks. The lag 24 result is fairly close to .05 and should not be described as strong proof of independence. The current TBATS script uses the default degrees of freedom in `Box.test`, so the reported p-values must be understood as diagnostics produced by that implementation.

Figure A3 gives additional evidence. One ACF lag was outside the bounds within the first 12 lags and two were outside within the first 24 lags. The residuals otherwise fluctuate around zero. Their distribution is close to a bell shape on the model's transformed scale. These isolated spikes show that the fit is acceptable but not perfect.

## III CONCLUSION

The fitted TBATS model provided the strongest overall forecasting accuracy. It ranked first in rolling-origin cross-validation and in the independent 2025 test set. Its small difference between training and CV errors also suggests stable generalisation. The residual tests passed at the .05 level. However, the lag 24 result and isolated ACF spikes require a cautious interpretation.

The model has two main limitations. Its automatically selected structure is less transparent than ETS or TSLM. It also relies only on past precipitation and cannot explain why an unusually wet month occurs. NASA POWER represents a spatial grid estimate rather than a local rain-gauge reading. Short local extremes may therefore be different.

Future work should confirm the residual test with degrees of freedom based on the fitted TBATS components. The model could also be compared with local station data and evaluated over more than one final test year. Used with these cautions, the forecasts can support rainfall preparedness work connected with Sustainable Development Goal 13 (United Nations, n.d.).

## REFERENCES

De Livera, A. M., Hyndman, R. J., & Snyder, R. D. (2011). Forecasting time series with complex seasonal patterns using exponential smoothing. *Journal of the American Statistical Association, 106*(496), 1513-1527. https://doi.org/10.1198/jasa.2011.tm09771

Hyndman, R. J. (2024). *forecast: Forecasting functions for time series and linear models* (Version 8.23.0) [Computer software]. https://pkg.robjhyndman.com/forecast/

Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and practice* (3rd ed.). OTexts. https://otexts.com/fpp3/

NASA Langley Research Center. (n.d.). *NASA POWER: Prediction of worldwide energy resources*. https://power.larc.nasa.gov/

United Nations. (n.d.). *Sustainable Development Goal 13: Climate action*. https://sdgs.un.org/goals/goal13

## APPENDIX A TABLES AND FIGURES

**Table A1. Accuracy of TBATS and the seasonal naive benchmark**

| Model | Evaluation data | RMSE | MAE | MAPE | MASE |
|---|---|---:|---:|---:|---:|
| TBATS | Training set | 2.449 | 1.878 | 38.844 | 0.730 |
| TBATS | 2025 test set | 2.335 | 1.984 | 35.702 | 0.771 |
| Seasonal naive | Training set | 3.319 | 2.574 | 53.406 | 1.000 |
| Seasonal naive | 2025 test set | 2.694 | 2.424 | 43.173 | 0.942 |

**Table A2. Cross-validation and residual diagnostic results for TBATS**

| Measure | TBATS result |
|---|---:|
| CV RMSE | 2.495 |
| CV MASE | 0.783 |
| SD of CV MASE | 0.189 |
| Complete CV folds | 27 |
| CV RMSE divided by training RMSE | 1.019 |
| Training to CV MASE difference | 6.84% |
| Ljung-Box p-value at lag 12 | .524 |
| Ljung-Box p-value at lag 24 | .089 |
| Significant ACF lags within first 12 | 1 |
| Significant ACF lags within first 24 | 2 |

**Figure A1. STL decomposition of Kuala Lumpur monthly mean precipitation from January 1981 to December 2025.**

Source image `output/plots/group_summary/stl_tbats_decomposition.png`

**Figure A2. TBATS forecasts for the 2025 test period with 80 percent and 95 percent prediction intervals and actual observations.**

Source image `output/plots/group_summary/fc_tbats.png`

**Figure A3. Time plot, autocorrelation function and distribution of residuals from the TBATS model.**

Source image `output/plots/group_summary/resid_tbats.png`

## APPENDIX B PROGRAMMING SOURCE CODE

Insert the complete contents of `scripts/06_ChiaZY_tbats.R`.
