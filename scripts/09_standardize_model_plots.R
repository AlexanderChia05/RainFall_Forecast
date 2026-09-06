# Regenerates aligned residual and forecast plots for all four models.
# Each model script is self-contained; sourcing them here provides one-click
# regeneration without requiring a pre-populated R session.

source("scripts/03_ChanYH_ets.R")
source("scripts/04_StephQF_arima.R")
source("scripts/05_HamGQ_tslm.R")
source("scripts/06_ChiaZY_tbats.R")

message("Aligned model plots regenerated in output/plots/group_summary")
