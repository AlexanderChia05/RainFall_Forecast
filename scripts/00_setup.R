user_library <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_library)) {
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("fpp3", "tseries", "zoo", "forecast", "Kendall")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, lib = .libPaths()[1])

library(dplyr)
library(purrr)
library(tidyr)
library(fpp3)
library(tseries)

acf_out_of_bounds <- function(resid, lag.max = 12) {
  r <- na.omit(resid)
  n <- length(r)
  ci <- 1.96 / sqrt(n)
  a <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}
