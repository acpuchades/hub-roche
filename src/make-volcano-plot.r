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

datos <- read_tsv(input_path) |>
  rename_with(normalize_names)

if ("or" %in% names(datos)) {
  datos$effect <- log2(datos$or)
  effect_label <- "$log_{2}(OR)$"
} else if ("beta" %in% names(datos)) {
  datos$effect <- datos$beta
  effect_label <- "$\\beta$"
} else {
  stop("No effect size column found")
}

datos |>
  mutate(
    p_value = do.call(coalesce, c(NA, across(any_of(c("permpvalue", "pvalue", "p"))))),
    log_p = pmin(log10(p_value), 1.2 * quantile(log10(p_value), 0.99, na.rm = TRUE)),
  ) |>
  slice_head(n = 1, by = c(chrom, pos, alt)) |>
  drop_na(p_value, effect) |>
  mutate(q_fdr = p.adjust(p, method = "fdr")) |>
  ggplot(aes(effect, -log_p, color = q_fdr <= fdr_threshold)) +
  geom_point() +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
  labs(
    x = TeX(effect_label), y = TeX("$-log_{10}(P)$"),
    color = TeX(str_glue("$FDR < {fdr_threshold}$")),
  )
ggsave(output_path, width = 12, height = 8, dpi = 300)
