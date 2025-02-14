library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(latex2exp)

source("src/utils.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

datos <- read_tsv(input_path) |> rename_with(normalize_names)

datos |>
  drop_na(p, or) |>
  mutate(q_fdr = p.adjust(p, method = "fdr")) |>
  ggplot(aes(log2(or), -log10(p), color = q_fdr < 0.05)) +
  geom_point() +
  coord_cartesian(xlim = c(-15, 15)) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    x = TeX("$log_{2}(OR)$"),
    y = TeX("$-log_{10}(P)$"),
    color = TeX("$FDR < 0.05$"),
  )
ggsave(output_path, width = 10, height = 10, dpi = 300)
