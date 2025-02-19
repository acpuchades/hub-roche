library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
pca_path <- args[1]
output_path <- args[2]

pca_data <- read_tsv(pca_path)

pca_data_l <- pca_data |>
  rename(IID = `#IID`) |>
  pivot_longer(starts_with("PC"), names_to = "PC") |>
  mutate(PC = as.integer(str_remove(PC, "PC")))

pca_data_l |>
  mutate(
    ax = PC %% 2,
    group = (PC - 1) %/% 2,
    group_label = str_glue("PC{group * 2 + 1} vs PC{group * 2 + 2}")
  ) |>
  select(-PC) |>
  pivot_wider(names_from = "ax", values_from = "value", names_prefix = "ax") |>
  ggplot(aes(ax0, ax1)) +
  geom_point() +
  facet_wrap(~group_label) +
  labs(x = NULL, y = NULL)

ggsave(output_path, width = 10, height = 6, dpi = 300)
