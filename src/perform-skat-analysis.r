library(cli)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

library(vcfR)
library(SKAT)

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- coalesce(args[1], "ALS")
skat_method <- coalesce(args[2], "SKATO")
maf_threshold <- coalesce(as.numeric(args[3]), 0.01)
vcf_path <- coalesce(args[4], "output/filtered-annotated-variants/cadd_20.vcf.gz")
psam_path <- coalesce(args[5], "output/variant-analysis/input/samples.psam")
cov_path <- coalesce(args[6], "output/variant-analysis/ALL_PCA.eigenvec")

cli_inform("Reading variants from {.file {vcf_path}}")
vcf_data <- read.vcfR(vcf_path, verbose = FALSE)
vcf_info <- getINFO(vcf_data)

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

sink(file = "/dev/null")

null_model <- SKAT_Null_Model(
  pheno_data[[pheno_name]] ~ PC1 + PC2 + PC3,
  data = pheno_data, out_type = "D", n.Resampling = 1000
)

results <- data.frame()
variant_sets <- str_extract(vcf_info, "SYMBOL=([^;]+)", group = 1)
for (name in unique(variant_sets)) {
  variant_idx <- which(variant_sets == name)
  set_variants <- gmat[variant_idx, , drop = FALSE]
  if (nrow(set_variants) == 1) {
    results <- rbind(results, data.frame(UNIT = name, P = NA))
    next
  }
  skat_model <- NULL
  try(skat_model <- SKATBinary(
    t(set_variants), null_model,
    method = skat_method,
    method.bin = "ER", max_maf = maf_threshold
  ))
  if (is.null(skat_model)) {
    results <- rbind(results, data.frame(UNIT = name, P = NA))
    next
  }
  results <- rbind(results, data.frame(
    UNIT = name, P = skat_model$p.value
  ))
}

sink()

write_tsv(results, stdout())
