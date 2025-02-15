library(readr)
library(dplyr)
library(HardyWeinberg)

source("src/utils.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

datos <- read_tsv(input_path) |>
  rename_with(normalize_names)

png(output_path, width = 800, height = 800)
HWTernaryPlot(
  as.matrix(datos[, 5:7]),
  vertexlab = c("AA", "AB", "BB")
)
dev.off()
