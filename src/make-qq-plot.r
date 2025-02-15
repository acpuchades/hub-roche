library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

source("src/utils.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

resultados <- read_tsv(input_path) |>
  rename_with(normalize_names) |>
  drop_na(p) |>
  mutate(
    observed = -log10(p),
    expected = -log10((rank(p) / (n() + 1))),
    q_fdr = p.adjust(p, method = "fdr")
  )

ggplot(resultados, aes(expected, observed, color = q_fdr <= 0.05)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text(
    aes(label = if_else(q_fdr <= 0.05, name, "")),
    nudge_y = -0.1
  ) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    x = TeX("Expected: $-log_{10}(P)$"),
    y = TeX("Observed: $-log_{10}(P)$"),
  ) +
  theme(legend.position = "none")
ggsave(output_path, width = 10, height = 9, dpi = 300)
