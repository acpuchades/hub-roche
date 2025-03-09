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
source("src/edmus.r")

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

ms_fclinica_data <- read_excel("data/ufem/FClinica.xlsx") |>
  rename_with(normalize_names) |>
  rename(id_edmus = fc_n_edmus) |>
  mutate(
    across(c(fc_nhc, fc_sap), as.integer),
    sexo = fc_sexo |> case_match(0 ~ "M", 1 ~ "F")
  )

ms_scores <- read_csv("data/ufem/datos-msss.csv") |>
  select(patient_id, msss = oGMSSS, armss = gARMSS)

ms_info_fclinica <- bind_rows(
  ms_fclinica_data |> mutate(nhc = fc_sap),
  ms_fclinica_data |> mutate(nhc = fc_nhc),
) |>
  select(nhc, fc_sap, fc_nhc, id_edmus, sexo) |>
  distinct()

ms_info_edmus <- edmus_personal |>
  select(-ms_onset) |>
  left_join(edmus_diagnosis, by = "patient_id") |>
  mutate(
    edad_inicio = floor((ms_onset - date_of_birth) / dyears(1)),
    fenotipo = disease_course |> case_match(
      "RR" ~ "RR",
      "SP-R" ~ "SP", "SP-NR" ~ "SP",
      "PP-R" ~ "PP", "PP-NR" ~ "PP"
    ),
    tiempo_hasta_iedss_3 = if_else(
      !irreversible_dss_3_unknown_date,
      round(
        (irreversible_dss_3 - ms_onset) / dmonths(1)
      ),
      NA_integer_
    ),
    tiempo_hasta_iedss_6 = if_else(
      !irreversible_dss_6_unknown_date,
      round((irreversible_dss_6 - ms_onset) / dmonths(1)),
      NA_integer_
    )
  ) |>
  left_join(edmus_arr, by = "patient_id")

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
    .keep = "unused"
  ) |>
  left_join(ms_info_fclinica, by = "nhc", na_matches = "never") |>
  left_join(ms_info_edmus, by = c(id_edmus = "local_identifier"), na_matches = "never")

# edmus_clinical |>
#  select(patient_id, date, edss = "edss_entered") |>
#  left_join(edmus_personal |> select(patient_id, date_of_birth, ms_onset), by = "patient_id") |>
#  filter(date |> between(ms_onset + dyears(2.5), ms_onset + dyears(3.5))) |>
#  slice_min(date, by = patient_id, n = 1, with_ties = FALSE) |>
#  semi_join(ms_samples_info, by = "patient_id") |>
#  transmute(
#    patient_id, dd = floor((date - ms_onset) / dyears(1)),
#    ageatedss = floor((date - date_of_birth) / dyears(1)),
#    edss
#  ) |>
#  write_csv("data/ufem/input-msss.csv")

ms_scores <- read_csv("data/ufem/datos-msss.csv") |>
  select(patient_id, msss = oGMSSS, armss = gARMSS)

ms_samples_info |>
  arrange(sample_id, grupo) |>
  left_join(ms_scores, by = "patient_id", na_matches = "never") |>
  transmute(
    SAMPLE = sample_id,
    NHC = nhc,
    SEX = sexo,
    PHENO = grupo |> case_match("CONTROL" ~ "CONTROL", "EM" ~ "MS"),
    MS_PHENO = fenotipo,
    MS_ONSET = edad_inicio,
    MS_MSSS = msss,
    MS_ARMSS = armss,
    MS_ARR_Y1 = arr_y1,
    MS_ARR_Y3 = round(arr_y3, 1),
    MS_IEDSS3 = tiempo_hasta_iedss_3,
    MS_IEDSS6 = tiempo_hasta_iedss_6,
  ) |>
  write_tsv(stdout())
