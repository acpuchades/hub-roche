library(tibble)
library(readr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(writexl)

library(DBI)
library(RSQLite)

source("src/utils.r")

normalize_sample_id <- function(x) {
  x |>
    str_replace("-?ADN[0-9]", "") |>
    str_replace("-ambBED", "") |>
    str_replace_all("-", "")
}

sample_ids <- read_lines("output/renamed-variants/samples.original.txt") |>
  as_tibble_col("sample_id")

ms_sample_ids <- sample_ids |>
  filter(str_starts(sample_id, "EM"))

als_sample_ids <- sample_ids |>
  filter(str_starts(sample_id, "ELA")) |>
  mutate(normalized_sample_id = normalize_sample_id(sample_id))

als_samples_bb <- read_excel("data/biobanco/Muestras ELA.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  filter(tipo_muestra %in% c("13-ADN", "11-Buffy coat"))

als_samples_ufela <- read_excel("data/ufela/biobanco-20240201.xlsx") |>
  rename(codigo_donacion = "...2") |>
  rename_with(normalize_names)

ufela_db <- dbConnect(RSQLite::SQLite(), "data/ufela/formulario_2023-11-15.sqlite")
ufela_pacientes <- dbGetQuery(ufela_db, "SELECT * FROM pacientes") |>
  rename_with(normalize_names)
dbDisconnect(ufela_db)

ms_samples_bb <- read_excel("data/biobanco/Muestras EM.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  mutate(nhc = na_if(nhc, "-"))

ms_samples_isa <- read_excel("data/ufem/Muestras-20240201.xlsx") |>
  rename(sample_id = id, nhc = "...3") |>
  rename_with(normalize_names) |>
  drop_na(sample_id)

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
    nhc = coalesce(nhc_biobanco, nhc_isa),
    diagnostico = coalesce(diagnostico_bb, diagnostico_isa),
    grupo = case_when(
      str_to_lower(diagnostico) |> str_detect("control") ~ "Control",
      !is.na(nhc) ~ "EM"
    )
  )

als_samples_info <- als_sample_ids |>
  left_join(
    als_samples_bb |> transmute(
      normalized_sample_id = normalize_sample_id(codigo_caso_noraybanks),
      nhc_noraybanks = nhc
    ),
    by = "normalized_sample_id", multiple = "first"
  ) |>
  left_join(
    als_samples_bb |> transmute(
      normalized_sample_id = normalize_sample_id(codigo_donacion_recepcion),
      nhc_donacion = nhc
    ),
    by = "normalized_sample_id", multiple = "first"
  ) |>
  left_join(
    als_samples_ufela |> transmute(
      normalized_sample_id = normalize_sample_id(caso_noraybanks),
      nhc_ufela_noraybanks = sap, diagnostico_noraybanks = dtco
    ),
    by = "normalized_sample_id", multiple = "first"
  ) |>
  left_join(
    als_samples_ufela |> transmute(
      normalized_sample_id = normalize_sample_id(codigo_donacion),
      nhc_ufela_donacion = sap, diagnostico_donacion = dtco
    ),
    by = "normalized_sample_id", multiple = "first"
  ) |>
  mutate(
    nhc = coalesce(nhc_noraybanks, nhc_donacion, nhc_ufela_noraybanks, nhc_ufela_donacion)
  ) |>
  left_join(
    ufela_pacientes |> transmute(nhc, diagnostico_ufela_db = "ELA"),
    by = "nhc"
  ) |>
  mutate(
    diagnostico = coalesce(diagnostico_noraybanks, diagnostico_donacion, diagnostico_ufela_db),
    grupo = if_else(
      str_to_lower(diagnostico) |> str_detect("control"),
      "Control", "ELA"
    ),
  )

als_samples_info |>
  select(sample_id, grupo) |>
  bind_rows(
    ms_samples_info |> select(sample_id, grupo)
  ) |>
  transmute(
    FID = sample_id,
    IID = sample_id,
    GROUP = case_match(
      grupo,
      "Control" ~ 0,
      "ELA" ~ 1,
      "EM" ~ 2,
      .default = -9
    ),
    ALS = if_else(grupo == "ELA", 1, 0, missing = -9),
    MS = if_else(grupo == "EM", 1, 0, missing = -9)
  ) |>
  write_tsv(stdout())
