library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(latex2exp)

source("src/utils.r")

fdr_threshold <- config::get("fdr_threshold")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

datos <- read_tsv(input_path) |> rename_with(normalize_names)

datos |>
  slice_head(n = 1, by = c(chrom, pos, alt)) |>
  drop_na(p, or) |>
  mutate(q_fdr = p.adjust(p, method = "fdr")) |>
  ggplot(aes(log2(or), -log10(p), color = q_fdr <= fdr_threshold)) +
  geom_point() +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
  labs(
    x = TeX("$log_{2}(OR)$"),
    y = TeX("$-log_{10}(P)$"),
    color = TeX(str_glue("$FDR < {fdr_threshold}$")),
  )
ggsave(output_path, width = 12, height = 8, dpi = 300)
