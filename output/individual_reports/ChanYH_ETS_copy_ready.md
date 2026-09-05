# Forecasting Monthly Mean Precipitation in Kuala Lumpur Using ETS

## I METHODOLOGY

### A Data preparation and initial evidence

This report analyses the monthly mean precipitation rate for Kuala Lumpur at 3.1390 N and 101.6869 E. The data were obtained from the NASA POWER project maintained by the NASA Langley Research Center (n.d.). PRECTOTCORR was measured in millimetres per day. The series covered January 1981 to December 2025 and contained 540 monthly observations. Values coded as -999 were treated as missing. Linear interpolation was included as a safeguard before the series was arranged in chronological order. The processed series contained no remaining missing values.

The final 12 months from January to December 2025 were reserved as an independent test set. The earlier 528 months were used for model fitting and cross-validation. This chronological split was necessary because a random split would allow later rainfall observations to influence the estimation of an earlier period.

STL decomposition was used to examine the pattern of the series. The seasonal strength was 0.469 while the trend strength was 0.181. Figure A1 shows a clear annual pattern and a weaker long-run movement. A Ljung-Box test on the raw series at lag 24 gave a test statistic of 233.12 with p < .001. The raw observations were therefore not white noise. These findings supported a model that updates the seasonal pattern while avoiding an unnecessary trend component.

### B ETS specification and justification

The selected model was ETS(A,N,A). The first A represents additive errors. N means that the trend state was omitted. The final A represents additive seasonality with a 12-month period. This form follows the state-space approach to exponential smoothing described by Hyndman et al. (2008).

The model estimates the level smoothing parameter and seasonal smoothing parameter by maximum likelihood. The no-trend setting was chosen because the STL trend strength was only 0.181. Additive seasonality was used because the rainfall cycle is expressed as changes in the original mm/day scale. No transformation was applied. The model was implemented with the fable package in R (O'Hara-Wild et al., 2024).

### C Forecast evaluation

Accuracy was assessed by RMSE, MAE, MAPE and MASE. Lower values indicate better forecasts. MASE was emphasised because it compares errors with the in-sample seasonal naive scale and is easier to compare across forecasting methods (Hyndman & Athanasopoulos, 2021). MAPE was retained for completeness but was interpreted with care because small rainfall values can increase percentage errors.

A seasonal naive model served as the benchmark. It forecasts each month using the observation from the same month one year earlier. Rolling-origin cross-validation started with 360 months of training data. The origin moved forward by six months and each fitted model forecasted the next 12 months. This produced 27 complete folds within the training period. The separate 2025 test set was not included in those folds.

## II DATA ANALYSIS RESULTS AND DISCUSSION

Table A1 shows that ETS achieved a training RMSE of 2.436 and a training MASE of 0.738. Its cross-validation RMSE was 2.566 and its cross-validation MASE was 0.811. The CV RMSE was 1.053 times the training RMSE. The MASE difference was 9.02 percent when measured against the CV value. Both changes were modest. The results do not show serious overfitting.

The seasonal naive model recorded a CV RMSE of 3.179 and a CV MASE of 0.985. ETS lowered CV MASE by about 17.6 percent. It also reduced CV RMSE by about 19.3 percent. The model therefore added forecasting value beyond repeating the previous year's monthly values.

On the independent 2025 test set, ETS obtained an RMSE of 2.363, an MAE of 2.077 and a MASE of 0.807. The seasonal naive test RMSE was 2.694 and its test MASE was 0.942. ETS reduced test RMSE by about 12.3 percent and test MASE by about 14.3 percent. Its test MAPE was 37.290 percent compared with 43.173 percent for the benchmark. Among the four main models, ETS ranked second on test RMSE and test MASE. TBATS was slightly more accurate.

Figure A2 compares the 2025 forecasts with the actual values. The forecast reproduces the broad annual movement but is smoother than several observed peaks and troughs. All 12 actual values appear within the plotted 95 percent prediction interval. This indicates useful interval coverage for the selected test year. It does not guarantee the same coverage in every future year.

Residual diagnostics were used to check whether the fitted model left a systematic time pattern. The Ljung-Box p-values were .522 at lag 12 and .122 at lag 24. Both exceed .05. The null hypothesis of no residual autocorrelation was therefore not rejected. This is a diagnostic pass rather than proof that the residuals are perfectly independent. Figure A3 shows no ACF spike outside the bounds within the first 12 lags and one isolated spike within 24 lags. The residuals remain centred near zero but the histogram is mildly right-skewed because some high-rainfall months were underpredicted.

## III CONCLUSION

ETS(A,N,A) was appropriate for a series with a visible annual cycle and weak trend strength. The model passed the residual checks used in this analysis. It also outperformed the seasonal naive benchmark in cross-validation and in the independent 2025 test set. The close training and CV errors suggest acceptable generalisation.

The model still has limitations. Additive seasonality assumes that seasonal movements are expressed on a roughly constant absolute scale. Figure A1 suggests that rainfall variability became wider in some later periods. The model also uses a NASA POWER grid estimate rather than a rain gauge at one station. Local rainfall extremes may therefore differ from the values analysed here. Weather drivers such as humidity and climate indices were not included.

Future work should compare additive ETS with a transformed or multiplicative alternative. Forecasts should also be checked against local station records. These improvements would be useful for rainfall monitoring and preparedness activities related to Sustainable Development Goal 13 on climate action (United Nations, n.d.).

## REFERENCES

Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and practice* (3rd ed.). OTexts. https://otexts.com/fpp3/

Hyndman, R. J., Koehler, A. B., Ord, J. K., & Snyder, R. D. (2008). *Forecasting with exponential smoothing: The state space approach*. Springer. https://doi.org/10.1007/978-3-540-71918-2

NASA Langley Research Center. (n.d.). *NASA POWER: Prediction of worldwide energy resources*. https://power.larc.nasa.gov/

O'Hara-Wild, M., Hyndman, R. J., & Wang, E. (2024). *fable: Forecasting models for tidy time series* (Version 0.3.4) [Computer software]. https://fable.tidyverts.org/

United Nations. (n.d.). *Sustainable Development Goal 13: Climate action*. https://sdgs.un.org/goals/goal13

## APPENDIX A TABLES AND FIGURES

**Table A1. Accuracy of ETS and the seasonal naive benchmark**

| Model | Evaluation data | RMSE | MAE | MAPE | MASE |
|---|---|---:|---:|---:|---:|
| ETS | Training set | 2.436 | 1.900 | 41.325 | 0.738 |
| ETS | 2025 test set | 2.363 | 2.077 | 37.290 | 0.807 |
| Seasonal naive | Training set | 3.319 | 2.574 | 53.406 | 1.000 |
| Seasonal naive | 2025 test set | 2.694 | 2.424 | 43.173 | 0.942 |

**Table A2. Cross-validation and residual diagnostic results for ETS**

| Measure | ETS result |
|---|---:|
| CV RMSE | 2.566 |
| CV MASE | 0.811 |
| SD of CV MASE | 0.188 |
| Complete CV folds | 27 |
| CV RMSE divided by training RMSE | 1.053 |
| Training to CV MASE difference | 9.02% |
| Ljung-Box p-value at lag 12 | .522 |
| Ljung-Box p-value at lag 24 | .122 |
| Significant ACF lags within first 12 | 0 |
| Significant ACF lags within first 24 | 1 |

**Figure A1. STL decomposition of Kuala Lumpur monthly mean precipitation from January 1981 to December 2025.**

Source image `output/plots/eda/stl_decomposition.png`

**Figure A2. ETS forecasts for the 2025 test period with 80 percent and 95 percent prediction intervals and actual observations.**

Source image `output/plots/group_summary/fc_ets.png`

**Figure A3. Time plot, autocorrelation function and distribution of residuals from the ETS model.**

Source image `output/plots/group_summary/resid_ets.png`

## APPENDIX B PROGRAMMING SOURCE CODE

Insert the complete contents of `scripts/03_ChanYH_ets.R`.
