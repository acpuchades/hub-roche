library(cli)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

library(vcfR)
library(SKAT)

pheno_categorical <- c(
  "SEX", "ALS", "MS", "ALS_MS",
  "ALS_UMN", "ALS_LMN", "ALS_RESTRICTED",
  "ALS_DFS_SP", "ALS_DFS_FP", "MS_SP", "MS_PP"
)

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- coalesce(args[1], "ALS")
skat_method <- coalesce(args[2], "SKATO")
maf_threshold <- coalesce(as.numeric(args[3]), 0.01)
vcf_path <- coalesce(args[4], "output/filtered-annotated-variants/impact_HIGH.vcf.gz")
vset_path <- coalesce(args[5], "output/variant-analysis/input/SYMBOL.vset")
psam_path <- coalesce(args[6], "output/variant-analysis/input/samples.psam")
cov_path <- coalesce(args[7], "output/variant-analysis/ALL_PCA.eigenvec")

cli_inform("Reading variants from {.file {vcf_path}}")
vcf_data <- read.vcfR(vcf_path, verbose = FALSE)
vcf_locus <- getFIX(vcf_data)[, c("CHROM", "POS", "REF", "ALT")]
vcf_ids <- str_c(vcf_locus[, 1], ":", vcf_locus[, 2], "_", vcf_locus[, 3], "_", vcf_locus[, 4])

cli_inform("Reading variant set information from {.file {vset_path}}")
vset_data <- read_tsv(vset_path, col_names = c("ID", "UNIT")) |>
  mutate(UNIT = UNIT |> na_if("."))

cli_inform("Reading phenotype information from {.file {psam_path}}")
psam_data <- read_tsv(psam_path) |>
  rename(IID = `#IID`) |>
  mutate(across(-IID, ~ .x |> case_match(1 ~ 0, 2 ~ 1)))

cli_inform("Reading additional covariates from {.file {cov_path}}")
cov_data <- read_tsv(cov_path) |>
  rename(IID = `#IID`)

pheno_data <- psam_data |>
  left_join(cov_data, by = "IID")

gmat <- vcf_data |>
  extract.gt(element = "GT")

gmat[gmat == "0/0"] <- "0"
gmat[gmat == "0/1" | gmat == "1/0"] <- "1"
gmat[gmat == "1/1"] <- "2"
gmat <- matrix(as.integer(gmat), nrow = nrow(gmat), ncol = ncol(gmat))

cli_inform("Starting {skat_method} analysis for {pheno_name} (MAF < {maf_threshold})")

sink(file = "/dev/null", type = "output")

null_model <- NULL
try(null_model <- SKAT_Null_Model(
  pheno_data[[pheno_name]] ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
  out_type = if_else(pheno_name %in% pheno_categorical, "D", "C"),
  data = pheno_data, n.Resampling = 1000
))

if (is.null(null_model)) {
  sink()
  write_tsv(data.frame(ID = NA, N = NA, P = NA), stdout())
  quit()
}

variant_sets <- tibble(ID = vcf_ids) |>
  left_join(vset_data, by = "ID")

results <- data.frame()
for (name in unique(variant_sets$UNIT)) {
  variant_idx <- which(variant_sets$UNIT == name)
  if (length(variant_idx) == 0) next
  set_variants <- gmat[variant_idx, , drop = FALSE]
  if (nrow(set_variants) == 1) {
    results <- rbind(results, data.frame(UNIT = name, N = 1, P = NA))
    next
  }
  skat_model <- NULL
  if (pheno_name %in% pheno_categorical) {
    try(skat_model <- SKATBinary(
      t(set_variants), null_model,
      method = skat_method, method.bin = "ER",
      max_maf = maf_threshold
    ))
  } else {
    try(skat_model <- SKAT(
      t(set_variants), null_model,
      method = skat_method, max_maf = maf_threshold
    ))
  }
  if (is.null(skat_model)) {
    results <- rbind(results, data.frame(UNIT = name, N = nrow(set_variants), P = NA))
    next
  }
  results <- rbind(results, data.frame(
    UNIT = name,
    N = nrow(set_variants),
    P = skat_model$p.value
  ))
}

sink(type = "output")
write_tsv(results, stdout())
