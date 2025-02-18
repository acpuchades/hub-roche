library(readr)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(latex2exp)

source("src/utils.r")

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
    log_p = -log10(p_value),
    q_fdr = p.adjust(p_value, method = "fdr")
  )

resultados |>
  slice_min(p_value, n = 20, with_ties = FALSE) |>
  mutate(name = fct_reorder(name, log_p)) |>
  ggplot(aes(x = log_p, y = name, fill = q_fdr <= config::get("fdr_threshold"))) +
  geom_col() +
  geom_text(
    aes(x = log_p / 2, label = round(q_fdr, digits = 3)),
    color = "white"
  ) +
  scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "red")) +
  labs(x = TeX("$-log_{10}(P)$"), y = NULL) +
  theme(legend.position = "none")
ggsave(output_path, width = 10, height = 5, dpi = 300)
