user_library <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_library)) {
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("fpp3", "tseries", "zoo", "forecast", "Kendall", "patchwork")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, lib = .libPaths()[1])

library(dplyr)
library(purrr)
library(tidyr)
library(fpp3)
library(tseries)
library(patchwork)

acf_out_of_bounds <- function(resid, lag.max = 12) {
  r <- na.omit(resid)
  n <- length(r)
  ci <- 1.96 / sqrt(n)
  a <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# Shared plot builders used by scripts 03-09.  Keeping these here ensures that
# every model script can run independently while producing the same layout.
residual_diagnostic_plot <- function(resid, index, model_name,
                                     acf_lag_max = 27) {
  d <- data.frame(t = index, r = as.numeric(resid))
  a <- acf(d$r, plot = FALSE, lag.max = acf_lag_max)$acf[-1]
  z <- data.frame(lag = seq_len(acf_lag_max), acf = a)
  ci <- 1.96 / sqrt(nrow(d))
  q <- theme_minimal() + theme(panel.grid.minor = element_blank())

  p1 <- ggplot(d, aes(t, r)) +
    geom_hline(yintercept = 0) + geom_line() + geom_point(size = 0.5) +
    labs(title = paste("Residual diagnostics -", model_name),
         x = NULL, y = "Residuals") + q
  p2 <- ggplot(z, aes(lag, acf)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-ci, ci), colour = "blue", linetype = "dashed") +
    geom_segment(aes(xend = lag, y = 0, yend = acf)) +
    scale_x_continuous(breaks = c(6, 12, 18, 24),
                       limits = c(1, acf_lag_max), expand = c(0, 0)) +
    labs(title = "ACF", x = "Lag [1M]", y = "ACF") + q
  p3 <- ggplot(d, aes(r)) +
    geom_histogram(aes(y = after_stat(density)), bins = 24, fill = "grey45") +
    stat_function(fun = dnorm,
                  args = list(mean = mean(d$r), sd = sd(d$r)),
                  colour = "orange") +
    labs(title = "Residual distribution", x = "Residuals", y = "Density") + q
  p1 / (p2 | p3)
}

save_residual_diagnostic <- function(resid, index, model_name, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, residual_diagnostic_plot(resid, index, model_name),
         width = 8, height = 6, dpi = 150)
}

save_fable_forecast_plot <- function(fc, actual, model_name, title, filename,
                                     plot_from = yearmonth("2013 Jan")) {
  test_actual <- actual |>
    as_tibble() |>
    filter(month > max(month) - 12) |>
    transmute(month, precip)
  p <- fc |>
    autoplot(actual |> filter(month >= plot_from), level = c(80, 95)) +
    geom_line(data = test_actual, aes(x = month, y = precip),
              colour = "red", linewidth = 0.45) +
    labs(title = title, y = "mm/day", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, p, width = 8, height = 5, dpi = 150)
  invisible(p)
}

save_tbats_forecast_plot <- function(fc_tbats, actual, train, filename,
                                     title = "TBATS: Forecast vs Actual",
                                     plot_from = yearmonth("2013 Jan")) {
  test_start <- max(train$month) + 1
  hist_tbl <- actual |> as_tibble() |>
    filter(month >= plot_from, month < test_start) |>
    transmute(month, value = precip)
  test_actual <- actual |> as_tibble() |>
    filter(month >= test_start) |>
    transmute(month, precip)
  h <- length(fc_tbats$mean)
  fc_tbl <- tibble(
    month = actual$month[(nrow(actual) - h + 1):nrow(actual)],
    value = as.numeric(fc_tbats$mean),
    lo80 = fc_tbats$lower[, "80%"], hi80 = fc_tbats$upper[, "80%"],
    lo95 = fc_tbats$lower[, "95%"], hi95 = fc_tbats$upper[, "95%"]
  )
  p <- ggplot() +
    geom_ribbon(data = fc_tbl, aes(month, ymin = lo95, ymax = hi95, fill = "95%"), alpha = .45) +
    geom_ribbon(data = fc_tbl, aes(month, ymin = lo80, ymax = hi80, fill = "80%"), alpha = .65) +
    geom_line(data = hist_tbl, aes(month, value), colour = "black", linewidth = .6) +
    geom_line(data = fc_tbl, aes(month, value), colour = "steelblue4", linewidth = .7) +
    geom_line(data = test_actual, aes(month, precip), colour = "red", linewidth = .45) +
    scale_fill_manual(name = "level", values = c("80%" = "#7F8CF0", "95%" = "#C5C8FF"),
                      breaks = c("80%", "95%")) +
    labs(title = title, y = "mm/day", x = NULL) +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, p, width = 8, height = 5, dpi = 150)
  invisible(p)
}
