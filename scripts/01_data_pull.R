source("scripts/00_setup.R")

url <- paste0(
  "https://power.larc.nasa.gov/api/temporal/monthly/point?",
  "parameters=PRECTOTCORR&community=AG",
  "&longitude=101.6869&latitude=3.1390",
  "&start=1981&end=2025&format=CSV"
)

raw_lines <- system(paste0("curl -s '", url, "'"), intern = TRUE)
start_line <- which(grepl("-END HEADER-", raw_lines)) + 1
df <- read.csv(text = paste(raw_lines[start_line:length(raw_lines)], collapse = "\n"),
                stringsAsFactors = FALSE)

month_abbr <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

rain <- df |>
  select(YEAR, all_of(month_abbr)) |>
  pivot_longer(-YEAR, names_to = "month_abbr", values_to = "precip") |>
  mutate(
    precip = na_if(precip, -999),
    month_num = match(month_abbr, month_abbr),
    month = yearmonth(paste(YEAR, month_num, sep = "-"))
  ) |>
  arrange(month) |>
  as_tsibble(index = month) |>
  mutate(precip = zoo::na.approx(precip, na.rm = FALSE)) |>
  select(month, precip)

dir.create("data", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
saveRDS(rain, "data/rain.rds")

cat("rain:", nrow(rain), "obs,", format(min(rain$month)), "to", format(max(rain$month)),
    "| NAs remaining:", sum(is.na(rain$precip)), "\n")
