source("scripts/00_setup.R")
rain <- readRDS("data/rain.rds")
dir.create("output/plots/eda", recursive = TRUE, showWarnings = FALSE)

# EDA
p_raw <- rain |> autoplot(precip) +
  labs(title = "KL monthly mean precipitation rate (mm/day)", y = "mm/day")
ggsave("output/plots/eda/raw_series.png", p_raw, width = 8, height = 5, dpi = 150)

p_box <- rain |>
  ggplot(aes(x = "Monthly precipitation", y = precip)) +
  geom_boxplot(fill = "steelblue", alpha = 0.65, width = 0.35) +
  labs(title = "Distribution of monthly precipitation", x = NULL, y = "mm/day") +
  theme_minimal()
ggsave("output/plots/eda/precipitation_boxplot.png", p_box, width = 6, height = 5, dpi = 150)

p_stl <- rain |>
  model(stl = STL(precip)) |>
  components() |>
  autoplot() +
  labs(title = "STL decomposition of monthly precipitation")
ggsave("output/plots/eda/stl_decomposition.png", p_stl, width = 8, height = 7, dpi = 150)

h <- 12
test_start <- max(rain$month) - h + 1
split_data <- rain |>
  mutate(data_set = factor(
    if_else(month < test_start, "Training set", "Test set"),
    levels = c("Training set", "Test set")
  ))

p_split <- ggplot(split_data, aes(month, precip, color = data_set)) +
  geom_line(linewidth = 0.55) +
  geom_vline(xintercept = as.Date(test_start), linetype = "dashed",
             color = "grey35") +
  scale_color_manual(values = c("Training set" = "#0072B2", "Test set" = "#D55E00")) +
  labs(title = "Training and test data", x = NULL, y = "mm/day", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("output/plots/eda/train_test_split.png", p_split, width = 8, height = 5, dpi = 150)

rain |> gg_season(precip) + labs(title = "Seasonal plot - by year")
rain |> gg_subseries(precip) + labs(title = "Subseries plot - by calendar month")
rain |> ACF(precip, lag_max = 36) |> autoplot() + labs(title = "ACF - raw series")
rain |> PACF(precip, lag_max = 36) |> autoplot() + labs(title = "PACF - raw series")

# Statistical tests
cat("\n== ADF test ==\n")
print(adf.test(rain$precip))

cat("\n== KPSS test ==\n")
print(kpss.test(rain$precip))

cat("\n== Ljung-Box test: raw series ==\n")
print(Box.test(rain$precip, lag = 12, type = "Ljung-Box"))
print(Box.test(rain$precip, lag = 24, type = "Ljung-Box"))

cat("\n== Mann-Kendall trend test ==\n")
print(Kendall::MannKendall(rain$precip))

cat("\n== Minimum precipitation ==\n")
print(min(rain$precip, na.rm = TRUE))
