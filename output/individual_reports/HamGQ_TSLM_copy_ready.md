# Forecasting Monthly Mean Precipitation in Kuala Lumpur Using TSLM

## I METHODOLOGY

### A Data preparation and exploratory findings

This study used monthly mean precipitation rates for Kuala Lumpur at 3.1390 N and 101.6869 E. The observations came from the NASA POWER database and were measured in millimetres per day (NASA Langley Research Center, n.d.). The series contained 540 months from January 1981 to December 2025. Values equal to -999 were defined as missing. Linear interpolation was included before the data were converted to a monthly time series. No missing observations remained after preparation.

The data were split chronologically. The model was fitted with 528 observations ending in December 2024. The 12 months in 2025 were held out for final testing. This approach preserved the order of the time series and prevented future values from entering the training process.

The STL analysis produced a seasonal strength of 0.469 and a trend strength of 0.181. Figure A1 shows that the repeating annual pattern is more visible than the gradual long-run movement. The raw series was also tested with ADF, KPSS, Mann-Kendall and Ljung-Box procedures. The Mann-Kendall test was used as a nonparametric supporting check for monotonic change (Mann, 1945). The raw Ljung-Box test at lag 24 gave a statistic of 233.12 with p < .001. The series therefore contained predictable time structure.

### B TSLM specification and justification

The selected model was `TSLM(precip ~ trend() + fourier(K = 4))`. It includes one linear time trend and four Fourier pairs for a 12-month seasonal cycle. Each Fourier pair contains one sine term and one cosine term. The model therefore used eight seasonal regressors together with the intercept and trend coefficient.

Fourier terms were selected because they describe a smooth repeating seasonal shape with fewer coefficients than 11 monthly indicators. K = 4 was used as a moderate level of seasonal detail. The trend term allowed the model to represent the weak long-run movement found in the exploratory analysis. All coefficients were estimated by ordinary least squares through the fable package (O'Hara-Wild et al., 2024). TSLM does not include an autoregressive error process in this specification. Residual checks were therefore important for deciding whether the trend and Fourier terms were sufficient.

### C Evaluation procedure

Forecasts were evaluated through training accuracy, rolling-origin cross-validation and a final test set. The CV procedure began with 360 observations. The origin advanced by six months and each model forecasted the next 12 months. A total of 27 complete folds were evaluated. The 2025 holdout remained separate.

RMSE, MAE, MAPE and MASE were reported. RMSE gives more weight to large forecasting mistakes. MAE describes the typical absolute error in mm/day. MASE compares the model with a seasonal naive error scale and is suitable for time series comparison (Hyndman & Athanasopoulos, 2021). The seasonal naive method was used as the benchmark.

## II DATA ANALYSIS RESULTS AND DISCUSSION

TSLM recorded a training RMSE of 2.417 and a training MASE of 0.736. The training RMSE was the lowest among the four main models. Training accuracy alone was not used to select the final model because a flexible regression can fit past values without producing the best future forecasts.

The mean CV RMSE was 2.557 and the mean CV MASE was 0.814. The CV RMSE was 1.058 times the training RMSE. The MASE difference was 9.56 percent relative to the CV value. This was the largest generalisation gap among the four main models but remained modest. The result suggests limited rather than serious overfitting.

The seasonal naive benchmark produced a CV RMSE of 3.179 and a CV MASE of 0.985. TSLM reduced CV RMSE by about 19.6 percent and CV MASE by about 17.3 percent. Its CV MASE was below one, so the regression provided useful improvement over repeating last year's monthly rainfall.

The independent 2025 test RMSE was 2.407. Test MAE was 2.146 and test MASE was 0.833. The test MAPE was 38.234 percent. All four errors were lower than the seasonal naive results. TSLM ranked third among the four main models on test RMSE and MASE. TBATS and ETS were more accurate. This difference shows why the lowest training RMSE should not be presented as evidence that TSLM was the best forecasting model.

Figure A2 shows that the fitted seasonal shape follows the general rise and fall of rainfall during 2025. The predicted line is smoother than the actual line because the seasonal component is formed from deterministic sine and cosine terms. Sharp monthly peaks are therefore difficult to reproduce. All 12 actual observations appear inside the plotted 95 percent prediction interval.

The Ljung-Box p-values were .476 at lag 12 and .148 at lag 24. Both values were above .05. The null hypothesis of no residual autocorrelation was not rejected. Among the four main models, TSLM had the largest lag 24 p-value. Figure A3 shows no significant ACF spike within the first 12 lags and one isolated spike within the first 24 lags. The residuals fluctuate around zero but the histogram has some right skew. The diagnostic evidence is acceptable even though the model has no explicit autoregressive error term.

## III CONCLUSION

The TSLM with a linear trend and four Fourier pairs was correctly implemented and produced forecasts that were better than the seasonal naive benchmark. Its residual checks passed at the .05 level. The difference between training and CV accuracy was small enough to indicate acceptable generalisation. However, it did not achieve the lowest CV or test errors.

The main strength of this model is its transparent structure. Trend and seasonal effects are represented directly through regression terms. Its limitation is that the seasonal shape is fixed and smooth. It may miss sudden monthly rainfall peaks or changes in monsoon timing. The model also contains no weather predictors and uses a NASA POWER grid estimate rather than a local station measurement.

Future research could add an ARIMA error component if residual dependence becomes stronger in new data. Climate variables could also be added as regressors. Validation with local rain-gauge records would improve practical confidence. These steps would make the analysis more useful for climate preparedness under Sustainable Development Goal 13 (United Nations, n.d.).

## REFERENCES

Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and practice* (3rd ed.). OTexts. https://otexts.com/fpp3/

Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica, 13*(3), 245-259. https://doi.org/10.2307/1907187

NASA Langley Research Center. (n.d.). *NASA POWER: Prediction of worldwide energy resources*. https://power.larc.nasa.gov/

O'Hara-Wild, M., Hyndman, R. J., & Wang, E. (2024). *fable: Forecasting models for tidy time series* (Version 0.3.4) [Computer software]. https://fable.tidyverts.org/

United Nations. (n.d.). *Sustainable Development Goal 13: Climate action*. https://sdgs.un.org/goals/goal13

## APPENDIX A TABLES AND FIGURES

**Table A1. Accuracy of TSLM with Fourier terms and the seasonal naive benchmark**

| Model | Evaluation data | RMSE | MAE | MAPE | MASE |
|---|---|---:|---:|---:|---:|
| TSLM with Fourier terms | Training set | 2.417 | 1.895 | 42.440 | 0.736 |
| TSLM with Fourier terms | 2025 test set | 2.407 | 2.146 | 38.234 | 0.833 |
| Seasonal naive | Training set | 3.319 | 2.574 | 53.406 | 1.000 |
| Seasonal naive | 2025 test set | 2.694 | 2.424 | 43.173 | 0.942 |

**Table A2. Cross-validation and residual diagnostic results for TSLM with Fourier terms**

| Measure | TSLM result |
|---|---:|
| CV RMSE | 2.557 |
| CV MASE | 0.814 |
| SD of CV MASE | 0.176 |
| Complete CV folds | 27 |
| CV RMSE divided by training RMSE | 1.058 |
| Training to CV MASE difference | 9.56% |
| Ljung-Box p-value at lag 12 | .476 |
| Ljung-Box p-value at lag 24 | .148 |
| Significant ACF lags within first 12 | 0 |
| Significant ACF lags within first 24 | 1 |

**Figure A1. STL decomposition of Kuala Lumpur monthly mean precipitation from January 1981 to December 2025.**

Source image `output/plots/eda/stl_decomposition.png`

**Figure A2. TSLM with Fourier forecasts for the 2025 test period with 80 percent and 95 percent prediction intervals and actual observations.**

Source image `output/plots/group_summary/fc_tslm.png`

**Figure A3. Time plot, autocorrelation function and distribution of residuals from the TSLM with Fourier model.**

Source image `output/plots/group_summary/resid_tslm.png`

## APPENDIX B PROGRAMMING SOURCE CODE

Insert the complete contents of `scripts/05_HamGQ_tslm.R`.
