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

sample_ids <- read_lines("output/renamed-variants/samples.normalized.txt") |>
  as_tibble_col("sample_id")

ms_sample_ids <- sample_ids |>
  filter(str_starts(sample_id, "EM"))

als_sample_ids <- sample_ids |>
  filter(str_starts(sample_id, "ELA"))

ms_samples_bb <- read_excel("data/biobanco/Muestras EM.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  mutate(nhc = na_if(nhc, "-"))

ms_samples_isa <- read_excel("data/ufem/Muestras-20240201.xlsx") |>
  rename(sample_id = id, nhc = "...3") |>
  rename_with(normalize_names) |>
  drop_na(sample_id)

ms_info_fclinica <- read_excel("data/ufem/FClinica.xlsx") |>
  rename_with(normalize_names) |>
  mutate(sexo = fc_sexo |> case_match(0 ~ "M", 1 ~ "F"))

ms_samples_info <- ms_sample_ids |>
  left_join(
    ms_samples_bb |> transmute(
      sample_id = adn,
      nhc_biobanco = nhc,
      diagnostico_bb = diagnostic
    ),
    by = "sample_id"
  ) |>
  left_join(
    ms_samples_isa |> transmute(
      sample_id,
      nhc_isa = nhc,
      diagnostico_isa = if_else(is.na(nhc), "CONTROL", "EM")
    ),
    by = "sample_id"
  ) |>
  mutate(
    nhc = as.integer(coalesce(nhc_biobanco, nhc_isa)),
    diagnostico = coalesce(diagnostico_bb, diagnostico_isa),
    grupo = case_when(
      str_to_lower(diagnostico) |> str_detect("control") ~ "Control",
      !is.na(nhc) ~ "EM"
    )
  ) |>
  left_join(ms_info_fclinica |> select(nhc = "fc_nhc", sexo), by = "nhc", na_matches = "never") |>
  left_join(ms_info_fclinica |> select(nhc = "fc_sap", sexo), by = "nhc", na_matches = "never") |>
  mutate(sexo = coalesce(sexo.x, sexo.y), .keep = "unused")

als_biobank_samples <- read_excel("data/biobanco/Muestras ELA.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  transmute(
    suppressWarnings(across(nhc, as.integer)),
    across(c(codigo_caso_noraybanks, codigo_donacion_recepcion), normalize_biobank_id),
    codigo_donacion_recepcion = biobank_donor_id(codigo_donacion_recepcion),
  ) |>
  distinct()

als_dropbox_samples <- suppressMessages(read_excel("data/ufela/biobanco-20240201.xlsx")) |>
  rename(codigo_donacion_recepcion = `...2`) |>
  rename_with(normalize_names) |>
  rename(codigo_caso_noraybanks = caso_noraybanks) |>
  mutate(nhc = suppressWarnings(as.integer(sap))) |>
  select(nhc, codigo_caso_noraybanks, codigo_donacion_recepcion) |>
  distinct()

als_ufela_db <- dbConnect(RSQLite::SQLite(), "data/ufela/formulario_2023-11-15.sqlite")
als_patients <- dbGetQuery(als_ufela_db, "SELECT * FROM pacientes") |>
  rename_with(normalize_names) |>
  mutate(nhc = suppressWarnings(as.integer(nhc)))
dbDisconnect(als_ufela_db)

als_ufela_samples <- als_biobank_samples |>
  bind_rows(als_dropbox_samples) |>
  distinct()

als_nhc_ccn <- als_ufela_samples |>
  select(nhc, codigo_caso_noraybanks) |>
  drop_na() |>
  distinct()

als_nhc_cdr <- als_ufela_samples |>
  select(nhc, codigo_donacion_recepcion) |>
  drop_na() |>
  distinct()

als_samples_nhc <- als_sample_ids |>
  left_join(
    bind_rows(
      als_nhc_ccn |> transmute(sample_id = codigo_caso_noraybanks, nhc, biobank_tag = "ccn"),
      als_nhc_cdr |> transmute(sample_id = codigo_donacion_recepcion, nhc, biobank_tag = "cdr")
    ),
    by = "sample_id", na_matches = "never"
  )

als_groups <- bind_rows(
  read_excel("data/ufela/grupos-20240112.xlsx", col_names = "sample_id", sheet = "ELA") |> mutate(grupo = "ELA"),
  read_excel("data/ufela/grupos-20240112.xlsx", col_names = "sample_id", sheet = "Controles") |> mutate(grupo = "Control"),
)

als_groups_nhc <- als_groups |>
  left_join(
    als_ufela_samples |> select(sample_id = codigo_donacion_recepcion, nhc),
    by = "sample_id", na_matches = "never"
  ) |>
  distinct()

als_samples_info <- als_samples_nhc |>
  left_join(als_groups_nhc |> select(nhc, grupo), by = "nhc", na_matches = "never") |>
  left_join(als_patients |> transmute(nhc, sexo, grupo = "ELA"), by = "nhc") |>
  mutate(
    grupo = coalesce(grupo.x, grupo.y), .keep = "unused",
    sexo = case_match(sexo, "Hombre" ~ "M", "Mujer" ~ "F")
  ) |>
  rows_update(tibble(sample_id = "ELA191", nhc = 18743408, sexo = "M", grupo = "ELA"), by = "sample_id") |>
  rows_update(tibble(sample_id = "ELA215", nhc = 17070408, sexo = NA, grupo = NA), by = "sample_id") |>
  rows_update(tibble(sample_id = "ELA218", nhc = NA, sexo = NA, grupo = NA), by = "sample_id") |>
  rows_update(tibble(sample_id = "ELA231", nhc = NA, sexo = NA, grupo = NA), by = "sample_id") |>
  rows_update(tibble(sample_id = "ELA233", nhc = 16844241, sexo = "F", grupo = "ELA"), by = "sample_id") |>
  rows_update(tibble(sample_id = "ELA236", nhc = 13783143, sexo = "F", grupo = "ELA"), by = "sample_id") |>
  slice_sample(n = 1, by = sample_id)

bind_rows(
  als_samples_info |> select(sample_id, grupo, sexo),
  ms_samples_info |> select(sample_id, grupo, sexo)
) |>
  arrange(sample_id) |>
  transmute(
    SAMPLE = sample_id,
    SEX = sexo,
    PHENO = grupo |> case_match("Control" ~ "CONTROL", "EM" ~ "MS", "ELA" ~ "ALS")
  ) |>
  write_tsv(stdout())
