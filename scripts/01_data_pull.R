# 01_data_pull.R - SHARED TOPIC (v3): monthly mean precipitation rate
# (mm/day) at Kuala Lumpur (3.139N, 101.6869E), NASA POWER (MERRA-2
# reanalysis). SDG 13 (Climate Action, indicator 13.1) primary, SDG 11
# (11.5, urban flood resilience) and SDG 6 (6.4, water resource
# planning) as extensions - same overarching SDG portfolio the group
# already committed to for the ozone (v1) topic, primary/extension
# swapped since rainfall maps onto climate action more directly than
# responsible consumption.
#
# Source: https://power.larc.nasa.gov/ (public API, no key required)
# Format verified before pulling: 7 header lines ending "-END HEADER-",
# then one row per YEAR with columns YEAR,JAN..DEC,ANN (mm/day). -999 is
# NASA POWER's missing-value placeholder (same convention family as the
# -9999 used in the ozone project's NASA Ozone Watch source).

source("scripts/00_setup.R")

url <- paste0(
  "https://power.larc.nasa.gov/api/temporal/monthly/point?",
  "parameters=PRECTOTCORR&community=AG",
  "&longitude=101.6869&latitude=3.1390",
  "&start=1981&end=2025&format=CSV"
)
# base R's readLines(url) hung indefinitely against this endpoint on
# Posit (verified: curl itself reaches it in <1s) - shell out to curl
# directly instead, same connection, no hang.
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
