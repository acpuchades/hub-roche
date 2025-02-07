library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(tidyr)

library(ggplot2)
library(ggvenn)
library(scales)
library(patchwork)
library(GGally)

as_chromosome <- function(x) {
  factor(x, levels = c(1:22, "X", "Y", "MT"))
}

variant.stats <- read_tsv("output/analysis-report/variant-stats.tsv") |>
  rename_with(\(x) x |> str_replace("^#?\\[[0-9]+\\]", "")) |>
  mutate(
    across(everything(), \(x) na_if(as.character(x), ".")),
    across(c(DP, QUAL, starts_with("gno_ac"), starts_with("gno_af")), as.numeric),
    across(CHROM, as_chromosome),
    locus = str_c(CHROM, ":", POS)
  )

variants_per_locus <- variant.stats |>
  slice_sample(n = 1, by = c(locus, "REF", "ALT"))

variant.scores <- variants_per_locus |>
  select(locus, REF, ALT, ends_with("_rankscore"), matches("SpliceAI_pred_DS_[AD][GL]")) |>
  pivot_longer(-c(locus, REF, ALT), names_to = "score", values_to = "value", values_transform = as.numeric) |>
  mutate(score = score |>
    str_replace_all("_converted", "") |>
    str_replace_all("_raw", "") |>
    str_replace_all("_rankscore$", ""))

variant.scores |>
  filter(!str_starts(score, "SpliceAI_")) |>
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~score, scales = "free") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Impact prediction scores", x = "Rank score", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("output/analysis-report/variants-predicted.png")

variant.scores |>
  filter(str_starts(score, "SpliceAI_")) |>
  mutate(score = score |>
    str_replace("SpliceAI_pred_DS_", "") |>
    case_match(
      "AG" ~ "Acceptor Gain",
      "AL" ~ "Acceptor Loss",
      "DG" ~ "Donor Gain",
      "DL" ~ "Donor Loss"
    )) |>
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~score) +
  scale_x_continuous(trans = "log10") +
  labs(title = "SpliceAI delta scores", x = "Delta score", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("output/analysis-report/variants-spliceai.png")

af_plot <- ggplot(variants_per_locus, aes(gno_af_all)) +
  geom_histogram() +
  scale_x_continuous(labels = label_percent())

af_plot_log10 <- ggplot(variants_per_locus, aes(gno_af_all)) +
  geom_histogram() +
  scale_x_continuous(name = "gno_af_all", labels = label_percent(), trans = "log10")

(af_plot | af_plot_log10) + plot_annotation(
  title = "Allele frequency of variants according to gnomAD genomes",
  theme = theme(plot.title = element_text(hjust = 0.5)),
)

ggsave("output/analysis-report/gnomad-af.png", width = 10, height = 5)

variants_per_locus |>
  transmute(
    locus,
    gnomAD = !is.na(gno_id),
    dbSNP = !is.na(dbsnp_id),
    ClinVar = !is.na(clinvar_id),
  ) |>
  ggvenn()
ggsave("output/analysis-report/variants-reported.png")

variants_per_locus |>
  mutate(CHROM = fct_rev(CHROM)) |>
  ggplot(aes(DP, CHROM)) +
  geom_boxplot() +
  scale_x_continuous(trans = "log10", labels = label_number()) +
  labs(title = "Variant coverage per chromosome")
ggsave("output/analysis-report/variants-coverage.png")

variants_per_locus |>
  mutate(CHROM = fct_rev(CHROM)) |>
  ggplot(aes(QUAL, CHROM)) +
  geom_boxplot() +
  scale_x_continuous(trans = "log10", labels = label_number()) +
  labs(title = "Variant call quality per chromosome")
ggsave("output/analysis-report/variants-quality.png")
