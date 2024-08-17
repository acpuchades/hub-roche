library(dplyr)
library(readr)
library(stringr)
library(forcats)

library(ggplot2)
library(ggvenn)
library(scales)

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

ggplot(variants_per_locus, aes(gno_af_all)) +
  geom_histogram() +
  scale_x_continuous(labels = label_percent()) +
  labs(title = "Variant frequency according to gnomAD v4.1")
ggsave("output/analysis-report/gnomad-af.png")

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
