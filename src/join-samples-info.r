library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
als_samples_info_path <- coalesce(args[1], "output/phenotype-info/ALS.tsv")
ms_samples_info_path <- coalesce(args[2], "output/phenotype-info/MS.tsv")

als_samples_data <- read_tsv(als_samples_info_path)
ms_samples_data <- read_tsv(ms_samples_info_path)

bind_rows(als_samples_data, ms_samples_data) |>
  write_tsv(stdout())
