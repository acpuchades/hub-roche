library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(latex2exp)

source("src/utils.r")

fdr_threshold <- config::get("fdr_threshold")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

resultados <- read_tsv(input_path) |>
  rename_with(normalize_names) |>
  mutate(
    name = do.call(coalesce, c(NA, across(any_of(c("id", "gene", "unit"))))),
    p_value = do.call(coalesce, c(NA, across(any_of(c("permpvalue", "pvalue", "p")))))
  ) |>
  slice_head(n = 1, by = name) |>
  drop_na(p_value) |>
  mutate(
    observed = -log10(p_value),
    expected = -log10((rank(p_value) / (n() + 1))),
    q_fdr = p.adjust(p_value, method = "fdr")
  )

ggplot(resultados, aes(expected, observed, color = q_fdr <= fdr_threshold)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text(
    aes(label = if_else(q_fdr <= fdr_threshold, as.character(name), "")),
    nudge_y = -0.1
  ) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    x = TeX("Expected: $-log_{10}(P)$"),
    y = TeX("Observed: $-log_{10}(P)$"),
  ) +
  theme(legend.position = "none")
ggsave(output_path, width = 10, height = 9, dpi = 300)
