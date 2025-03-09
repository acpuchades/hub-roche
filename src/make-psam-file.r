library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
samples_info_path <- coalesce(args[1], "output/phenotype-info/ALL.tsv")

samples_info <- read_tsv(samples_info_path) |>
  arrange(SAMPLE) |>
  rename(`#IID` = SAMPLE) |>
  mutate(
    SEX = SEX |> case_match("M" ~ 1, "F" ~ 2),
    ALS = PHENO |> case_match("CONTROL" ~ 1, "ALS" ~ 2),
    MS = PHENO |> case_match("CONTROL" ~ 1, "MS" ~ 2),
    ALS_MS = PHENO |> case_match("CONTROL" ~ 1, "ALS" ~ 2, "MS" ~ 2),
    .after = SEX, .keep = "unused"
  ) |>
  mutate(
    ALS_UMN = ALS_PHENO |> case_match(
      "ALS-S" ~ 2, "ALS-B" ~ 2, "ALS-R" ~ 2,
      "PLS" ~ 2, "PMA" ~ 1,
    ),
    ALS_LMN = ALS_PHENO |> case_match(
      "ALS-S" ~ 2, "ALS-B" ~ 2, "ALS-R" ~ 2,
      "PLS" ~ 1, "PMA" ~ 2,
    ),
    ALS_RESTRICTED = ALS_PHENO |> case_match(
      "ALS-S" ~ 1, "ALS-B" ~ 1, "ALS-R" ~ 1,
      "PBP" ~ 2, "FAS" ~ 2, "FLS" ~ 2,
    ),
    ALS_DFS_SP = ALS_DFS_PC |> case_match("SP" ~ 2, "NP" ~ 1, "FP" ~ 1),
    ALS_DFS_FP = ALS_DFS_PC |> case_match("SP" ~ 1, "NP" ~ 1, "FP" ~ 2),
    ALS_D50_SP = ALS_D50_PC |> case_match("SP" ~ 2, "NP" ~ 1, "FP" ~ 1),
    ALS_D50_FP = ALS_D50_PC |> case_match("SP" ~ 1, "NP" ~ 1, "FP" ~ 2),
    .after = ALS_D50_PC, .keep = "unused"
  ) |>
  mutate(
    MS_SP = MS_PHENO |> case_match("RR" ~ 1, "SP" ~ 2, "PP" ~ 1),
    MS_PP = MS_PHENO |> case_match("RR" ~ 1, "SP" ~ 1, "PP" ~ 2),
    .after = MS_PHENO, .keep = "unused"
  ) |>
  select(-NHC)

write_tsv(samples_info, stdout())
