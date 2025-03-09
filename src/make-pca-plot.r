library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
pca_path <- coalesce(args[1], "output/variant-analysis/ALL_PCA.eigenvec")
samples_info_path <- coalesce(args[2], "output/phenotype-info/ALL.tsv")
output_path <- coalesce(args[3], "output/variant-analysis/ALL_PCA.eigenvec.png")

pca_data <- read_tsv(pca_path) |>
  rename(IID = `#IID`)

samples_info <- read_tsv(samples_info_path) |>
  select(IID = "SAMPLE", PHENO)

pca_data_l <- pca_data |>
  inner_join(samples_info |> select(IID, PHENO), by = "IID") |>
  pivot_longer(starts_with("PC"), names_to = "PC") |>
  mutate(PC = as.integer(str_remove(PC, "PC")) - 1)

pca_data_l |>
  mutate(group = PC %/% 2, ax = PC %% 2) |>
  select(-PC) |>
  pivot_wider(names_from = ax, names_prefix = "ax", values_from = value) |>
  mutate(group_label = str_glue("PC{group * 2 + 1} / PC{group * 2 + 2}")) |>
  ggplot(aes(ax0, ax1, color = PHENO)) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = NULL, y = NULL, color = NULL) +
  facet_wrap(~group_label) +
  scale_color_manual(values = c(
    `NA` = "black",
    ALS = "#F20403",
    MS = "#0B548B",
    CONTROL = "#138045"
  ))

ggsave(output_path, width = 10, height = 6, dpi = 300)
