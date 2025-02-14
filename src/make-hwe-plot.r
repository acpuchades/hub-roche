library(dplyr)
library(readr)

library(HardyWeinberg)

source("src/utils.r")

data <- read_tsv("output/variant-analysis/controls.hardy")

HWTernaryPlot(data[, 5:7],
  region = 0, markercol = rgb(0, 0, 1, 0.03),
)
