# 00_setup.R — packages. Run once per Posit Cloud session.
# Free tier 1GB RAM. Do NOT install prophet (Stan compile fails/times out).

pkgs <- c("fpp3", "tseries", "zoo", "forecast")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(dplyr)
library(purrr)
library(tidyr)
library(fpp3)
library(tseries)

# NOTE: forecast package is installed above (for TBATS/NNAR supplementary
# tests, 10_tbats_nnar_grid.R) but deliberately NOT attached with
# library() here - forecast::accuracy()/forecast() would mask
# fabletools::accuracy()/forecast(), which every other script in this
# project calls unqualified assuming the fabletools version. Scripts
# that need TBATS call forecast::tbats(), forecast::forecast(),
# forecast::accuracy() explicitly instead (same convention as the ozone
# project's BATS member).

# Shared topic (all 4 members): KL monthly mean precipitation (mm/day),
# NASA POWER (SDG 13 - see 01_data_pull.R header for full SDG writeup).
# acf_out_of_bounds(): count residual ACF lags (1..lag.max) exceeding the
# +-1.96/sqrt(n) critical value. Must-check per model - a handful (~5% of
# lag.max by chance) is fine, many = residual autocorrelation left behind.
#
# lag.max = 12 (one seasonal period, m), matching the convention carried
# over from the ozone project (V1) - standard alternative to 2m per
# Hyndman & Athanasopoulos. Verify both lag=12 and lag=24 per model
# rather than assuming this convention transfers cleanly to a new topic.
acf_out_of_bounds <- function(resid, lag.max = 12) {
  r <- na.omit(resid)
  n <- length(r)
  ci <- 1.96 / sqrt(n)
  a <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}
