# Forecasting Monthly Mean Precipitation in Kuala Lumpur Using ARIMA with Fourier Terms

## I METHODOLOGY

### A Data preparation and time series evidence

This report examines the monthly mean precipitation rate for Kuala Lumpur at 3.1390 N and 101.6869 E. The PRECTOTCORR data were obtained from NASA POWER and were measured in millimetres per day (NASA Langley Research Center, n.d.). The record covered January 1981 to December 2025 and contained 540 monthly observations. The preprocessing code converted -999 values to missing values and applied linear interpolation where required. The final ordered time series had no remaining missing observations.

The first 528 months ending in December 2024 formed the training set. January to December 2025 formed the 12-month test set. A chronological split was used so that future observations did not influence estimation or model selection.

STL decomposition gave a seasonal strength of 0.469 and a trend strength of 0.181. Figure A1 shows a repeating annual cycle with weaker long-run movement. The raw series also rejected the Ljung-Box white-noise test at lag 24 with a statistic of 233.12 and p < .001. This means that the observations contained time dependence that could potentially be forecasted.

### B ARIMA with Fourier specification

The selected specification was `ARIMA(precip ~ fourier(K = 4) + pdq() + PDQ(0,0,0))`. Four Fourier pairs represented the 12-month seasonal pattern. Each pair contains a sine and cosine term, giving eight deterministic seasonal regressors. K = 4 provides more detail than a simple single wave while avoiding a separate parameter for every month.

The seasonal ARIMA order was fixed at (0,0,0). This prevented seasonal ARIMA terms from competing with the Fourier regressors for the same annual pattern. The non-seasonal values of p, d and q were chosen automatically by the fable ARIMA search. The final error structure used one regular difference and one autoregressive parameter. It can therefore be written as ARIMA(1,1,0) errors with Fourier K = 4. The differencing diagnostic is shown in Figure A2. This approach follows the automatic ARIMA framework of Hyndman and Khandakar (2008). Model fitting used the fable package (O'Hara-Wild et al., 2024).

### C Evaluation procedure

The model was assessed by training accuracy, rolling-origin cross-validation and the independent 2025 test set. Cross-validation began with the first 360 months. The origin advanced by six months. Each model forecasted 12 months ahead. There were 27 complete folds and none used the 2025 observations.

RMSE, MAE, MAPE and MASE were reported. RMSE penalises large errors more strongly. MAE gives the average error size in mm/day. MASE compares the error with the seasonal naive scale and supports comparison across methods (Hyndman & Athanasopoulos, 2021). MAPE was treated as supporting evidence because percentage errors become less stable when actual rainfall is low.

## II DATA ANALYSIS RESULTS AND DISCUSSION

The ARIMA with Fourier model gave a training RMSE of 2.438 and a training MASE of 0.738. Its CV RMSE was 2.538 and its CV MASE was 0.801. The CV RMSE was 1.041 times the training RMSE. The MASE difference was 7.82 percent relative to the CV value. These small changes suggest that the model generalised well and did not show serious overfitting.

The seasonal naive model recorded a CV RMSE of 3.179 and a CV MASE of 0.985. ARIMA with Fourier terms reduced CV RMSE by about 20.2 percent and CV MASE by about 18.7 percent. Its CV MASE was the second lowest among the four main models. TBATS was slightly better with a CV MASE of 0.783.

On the independent 2025 test set, ARIMA obtained an RMSE of 2.424, an MAE of 2.170, a MAPE of 39.412 percent and a MASE of 0.843. The seasonal naive test RMSE was 2.694 and its test MASE was 0.942. ARIMA therefore improved test RMSE by about 10.0 percent and test MASE by about 10.5 percent. It passed the benchmark but ranked behind TBATS, ETS and TSLM on test RMSE and MASE. This does not make ARIMA unsuitable. It means that the latest evidence does not support calling it the overall best model.

Figure A2 displays the ACF and PACF after one regular difference. The remaining spikes explain why an ARIMA error structure was useful in addition to the deterministic seasonal terms. Figure A3 shows the final forecasts. The forecast follows the main annual movement but smooths several sharp changes in the actual 2025 series. All 12 observations appear inside the plotted 95 percent prediction interval.

The residual Ljung-Box test accounted for the one estimated AR or MA parameter. The p-value was .467 at lag 12 and .108 at lag 24. Both exceed .05. The null hypothesis of no remaining residual autocorrelation was not rejected. Figure A4 shows no significant ACF spike within the first 12 lags and one isolated spike within the first 24 lags. The histogram is centred near zero with a right tail caused by several large positive forecast errors. The diagnostic results support an adequate model but do not prove complete independence or normality.

## III CONCLUSION

ARIMA with four Fourier pairs successfully represented the annual pattern and the remaining short-run dependence. Its training and CV errors were close. It passed the residual diagnostics and performed better than the seasonal naive benchmark in both CV and the 2025 test set. It was a strong model but was not the overall winner. TBATS had lower CV and test errors.

The fixed Fourier terms assume that the seasonal shape remains stable. This may be unrealistic if monsoon timing changes. The model also smooths extreme monthly peaks. It uses only historical precipitation and cannot explain rainfall changes through physical weather factors. NASA POWER is a grid-based estimate rather than a direct local station record.

Future work could compare different K values through the same rolling-origin procedure instead of selecting K from training fit. Climate variables could be added through dynamic regression. Local station validation would also strengthen the result. These extensions would improve the model's value for climate preparedness under Sustainable Development Goal 13 (United Nations, n.d.).

## REFERENCES

Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and practice* (3rd ed.). OTexts. https://otexts.com/fpp3/

Hyndman, R. J., & Khandakar, Y. (2008). Automatic time series forecasting: The forecast package for R. *Journal of Statistical Software, 27*(3), 1-22. https://doi.org/10.18637/jss.v027.i03

NASA Langley Research Center. (n.d.). *NASA POWER: Prediction of worldwide energy resources*. https://power.larc.nasa.gov/

O'Hara-Wild, M., Hyndman, R. J., & Wang, E. (2024). *fable: Forecasting models for tidy time series* (Version 0.3.4) [Computer software]. https://fable.tidyverts.org/

United Nations. (n.d.). *Sustainable Development Goal 13: Climate action*. https://sdgs.un.org/goals/goal13

## APPENDIX A TABLES AND FIGURES

**Table A1. Accuracy of ARIMA with Fourier terms and the seasonal naive benchmark**

| Model | Evaluation data | RMSE | MAE | MAPE | MASE |
|---|---|---:|---:|---:|---:|
| ARIMA with Fourier terms | Training set | 2.438 | 1.901 | 41.193 | 0.738 |
| ARIMA with Fourier terms | 2025 test set | 2.424 | 2.170 | 39.412 | 0.843 |
| Seasonal naive | Training set | 3.319 | 2.574 | 53.406 | 1.000 |
| Seasonal naive | 2025 test set | 2.694 | 2.424 | 43.173 | 0.942 |

**Table A2. Cross-validation and residual diagnostic results for ARIMA with Fourier terms**

| Measure | ARIMA result |
|---|---:|
| CV RMSE | 2.538 |
| CV MASE | 0.801 |
| SD of CV MASE | 0.194 |
| Complete CV folds | 27 |
| CV RMSE divided by training RMSE | 1.041 |
| Training to CV MASE difference | 7.82% |
| Ljung-Box p-value at lag 12 | .467 |
| Ljung-Box p-value at lag 24 | .108 |
| Significant ACF lags within first 12 | 0 |
| Significant ACF lags within first 24 | 1 |

**Figure A1. STL decomposition of Kuala Lumpur monthly mean precipitation from January 1981 to December 2025.**

Source image `output/plots/eda/stl_decomposition.png`

**Figure A2. ACF and PACF of the Kuala Lumpur precipitation series after one regular difference.**

Source image `output/plots/group_summary/arima_differenced_acf_pacf.png`

**Figure A3. ARIMA with Fourier forecasts for the 2025 test period with 80 percent and 95 percent prediction intervals and actual observations.**

Source image `output/plots/group_summary/fc_arima.png`

**Figure A4. Time plot, autocorrelation function and distribution of residuals from the ARIMA with Fourier model.**

Source image `output/plots/group_summary/resid_arima.png`

## APPENDIX B PROGRAMMING SOURCE CODE

Insert the complete contents of `scripts/04_StephQF_arima.R`.
