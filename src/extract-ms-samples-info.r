library(tibble)
library(readr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(lubridate)
library(writexl)

library(DBI)
library(RSQLite)

source("src/utils.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- coalesce(args[1], "output/renamed-variant-samples/samples.normalized.txt")
vcf_samples <- read_lines(input_path) |> as_tibble_col("sample_id")

ms_sample_ids <- vcf_samples |>
  filter(str_starts(sample_id, "EM"))

ms_samples_bb <- read_excel("data/biobanco/Muestras EM.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  select(nhc, sample_id = "codi_bb", grupo = "diagnostic") |>
  drop_na(nhc, sample_id) |>
  mutate(
    across(nhc, as.integer),
    grupo = coalesce(grupo, "EM")
  )

ms_samples_isa <- read_excel("data/ufem/Muestras-20240201.xlsx") |>
  rename(sample_id = id, nhc = "...3") |>
  rename_with(normalize_names) |>
  drop_na(sample_id) |>
  mutate(
    across(nhc, as.integer),
    grupo = if_else(is.na(nhc), "CONTROL", "EM")
  )

ms_info_fclinica <- read_excel("data/ufem/FClinica.xlsx") |>
  rename_with(normalize_names) |>
  mutate(
    across(c(fc_nhc, fc_sap), as.integer),
    sexo = fc_sexo |> case_match(0 ~ "M", 1 ~ "F")
  )

ms_samples_info <- ms_sample_ids |>
  left_join(
    ms_samples_bb |> select(sample_id, nhc_bb = "nhc", grupo_bb = "grupo"),
    by = "sample_id"
  ) |>
  left_join(
    ms_samples_isa |> select(sample_id, nhc_isa = "nhc", grupo_isa = "grupo"),
    by = "sample_id"
  ) |>
  mutate(
    nhc = coalesce(nhc_bb, nhc_isa),
    grupo = coalesce(grupo_bb, grupo_isa),
  ) |>
  left_join(ms_info_fclinica |> select(nhc = "fc_nhc", sexo), by = "nhc", na_matches = "never") |>
  left_join(ms_info_fclinica |> select(nhc = "fc_sap", sexo), by = "nhc", na_matches = "never") |>
  mutate(sexo = coalesce(sexo.x, sexo.y), .keep = "unused")

ms_samples_info |>
  arrange(sample_id) |>
  transmute(
    SAMPLE = sample_id,
    NHC = nhc,
    SEX = sexo,
    PHENO = grupo |> case_match("CONTROL" ~ "CONTROL", "EM" ~ "MS")
  ) |>
  write_tsv(stdout())
