library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(latex2exp)

source("src/utils.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

datos <- read_tsv(input_path) |>
  rename_with(normalize_names) |>
  mutate(chrom_num = coalesce(
    as.integer(chrom),
    chrom |> case_match("X" ~ 23, "Y" ~ 24, "MT" ~ 25)
  ))

datos_acc <- datos |>
  group_by(chrom_num) |>
  summarize(max_pos = max(pos, na.rm = TRUE)) |>
  mutate(pos_add = lag(cumsum(max_pos), default = 0))

datos <- datos |>
  left_join(datos_acc, by = "chrom_num") |>
  mutate(pos_acc = pos + pos_add)

axis_set <- datos |>
  group_by(chrom) |>
  summarize(center = mean(pos_acc))

datos |>
  ggplot(aes(pos_acc, -log10(p), color = as.factor(chrom_num %% 2))) +
  geom_point(size = 1) +
  scale_x_continuous(
    label = axis_set$chrom,
    breaks = axis_set$center
  ) +
  scale_color_manual(values = c("#276FBF", "#183059")) +
  labs(x = "Chromosome", y = TeX("$-log_{10}(P)$"), color = NULL) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5)
  )
ggsave(output_path, width = 10, height = 5, dpi = 300)
