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
source("src/ufela.r")

args <- commandArgs(trailingOnly = TRUE)
input_path <- coalesce(args[1], "output/renamed-variant-samples/samples.normalized.txt")
sample_ids <- read_lines(input_path) |> as_tibble_col("sample_id")

als_sample_ids <- sample_ids |>
  filter(str_starts(sample_id, "ELA"))

als_samples_bb <- read_excel("data/biobanco/Muestras ELA.xlsx", skip = 1) |>
  rename_with(normalize_names) |>
  transmute(
    suppressWarnings(across(nhc, as.integer)),
    across(c(codigo_caso_noraybanks, codigo_donacion_recepcion), normalize_biobank_id),
    codigo_donante = biobank_donor_id(codigo_donacion_recepcion),
    fecha_muestra = dmy(coalesce(fecha_muestra, fecha_donacion_recepcion))
  )

als_dropbox_data <- suppressMessages(read_excel("data/ufela/biobanco-20240201.xlsx")) |>
  rename(codigo_donacion_recepcion = `...2`) |>
  rename_with(normalize_names) |>
  filter(codigo_donacion_recepcion != "Biobanc") |>
  rename(codigo_caso_noraybanks = caso_noraybanks) |>
  mutate(
    across(starts_with("fecha_"), ~ suppressWarnings(coalesce(
      dmy(.x), as_date(as.numeric(.x), origin = "1899-12-30")
    ))),
    nhc = suppressWarnings(as.integer(sap)),
    grupo = case_when(
      str_to_lower(dtco) |> str_detect("control") ~ "CONTROL",
      str_to_lower(dtco) |> str_detect("ela|elp|amp|pma|pbp|(flail arm|leg)") ~ "ELA",
    )
  )

als_samples_dropbox <- als_dropbox_data |>
  select(nhc, codigo_caso_noraybanks, codigo_donacion_recepcion, fecha_muestra)

als_samples <- als_samples_bb |>
  bind_rows(als_samples_dropbox) |>
  separate_longer_delim(cols = codigo_donacion_recepcion, delim = "/") |>
  mutate(
    across(c(codigo_caso_noraybanks, codigo_donacion_recepcion), normalize_biobank_id),
    codigo_donante = biobank_donor_id(codigo_donacion_recepcion)
  ) |>
  distinct()

als_groups_roche <- bind_rows(
  read_excel("data/ufela/grupos-20240112.xlsx", col_names = "codigo_donacion_recepcion", sheet = "ELA") |> mutate(grupo = "ELA"),
  read_excel("data/ufela/grupos-20240112.xlsx", col_names = "codigo_donacion_recepcion", sheet = "Controles") |> mutate(grupo = "Control"),
) |> mutate(grupo = grupo |> case_match("Control" ~ "CONTROL", "ELA" ~ "ELA"))

als_groups_dropbox <- als_dropbox_data |>
  select(codigo_caso_noraybanks, codigo_donacion_recepcion, grupo)

als_groups <- als_groups_roche |>
  bind_rows(als_groups_dropbox) |>
  mutate(codigo_donante = biobank_donor_id(codigo_donacion_recepcion))

als_ccn_nhc <- als_samples |>
  select(nhc, codigo_caso_noraybanks) |>
  drop_na() |>
  distinct() |>
  filter(n() == 1, .by = "codigo_caso_noraybanks")

als_cdr_nhc <- als_samples |>
  select(nhc, codigo_donante) |>
  drop_na() |>
  distinct() |>
  filter(n() == 1, .by = "codigo_donante")

als_samples_nhc <- als_sample_ids |>
  left_join(
    bind_rows(
      als_ccn_nhc |> transmute(sample_id = codigo_caso_noraybanks, nhc, biobank_tag = "ccn"),
      als_cdr_nhc |> transmute(sample_id = codigo_donante, nhc, biobank_tag = "cdr")
    ),
    by = "sample_id", na_matches = "never"
  )

als_ccn_group <- als_groups |>
  select(codigo_caso_noraybanks, grupo) |>
  drop_na() |>
  distinct() |>
  filter(n() == 1, .by = "codigo_caso_noraybanks")

als_cdr_group <- als_groups |>
  select(codigo_donante, grupo) |>
  drop_na() |>
  distinct() |>
  filter(n() == 1, .by = "codigo_donante")

als_nhc_group <- bind_rows(
  als_ccn_group |> inner_join(als_ccn_nhc, by = "codigo_caso_noraybanks", na_matches = "never"),
  als_cdr_group |> inner_join(als_cdr_nhc, by = "codigo_donante", na_matches = "never")
) |>
  pivot_longer(
    c(codigo_caso_noraybanks, codigo_donante),
    names_to = "biobank_tag", values_to = "sample_id"
  ) |>
  select(nhc, grupo) |>
  drop_na() |>
  distinct()

als_samples_info <- als_samples_nhc |>
  left_join(als_nhc_group, by = "nhc", na_matches = "never") |>
  # left_join(als_cdr_group |> select(sample_id = "codigo_donante", grupo), by = "sample_id", na_matches = "never") |>
  # mutate(grupo = coalesce(grupo.x, grupo.y), .keep = "unused") |>
  # left_join(als_ccn_group |> select(sample_id = "codigo_caso_noraybanks", grupo), by = "sample_id", na_matches = "never") |>
  # mutate(grupo = coalesce(grupo.x, grupo.y), .keep = "unused") |>
  left_join(ufela_datos |> transmute(
    nhc, sexo,
    grupo = "ELA",
    edad_inicio,
    retraso_diagnostico_meses,
    fenotipo = fenotipo_al_diagnostico |> case_match(
      "ELA Bulbar" ~ "ALS-B",
      "ELA Espinal" ~ "ALS-S",
      "Monomiélica" ~ "ALS-S",
      "ELA Respiratoria" ~ "ALS-R",
      "Parálisis bulbar progresiva" ~ "PBP",
      "Flail leg" ~ "FLS",
      "Flail arm" ~ "FAS",
      "Atrofia Muscular Progresiva (AMP)" ~ "PMA",
      "Esclerosis Lateral Primaria (ELP)" ~ "PLS",
    ),
    delta_fs = if_else(fecha_deltafs - fecha_diagnostico_ela <= dmonths(12), delta_fs, NA_real_),
    fvc_basal = if_else(fecha_fvc - fecha_diagnostico_ela <= dmonths(12), fvc_basal, NA),
    categoria_progresion = if_else(fecha_deltafs - fecha_diagnostico_ela <= dmonths(12), categoria_progresion, NA),
    supervivencia_meses = round((fecha_exitus - fecha_inicio_clinica) / dmonths(1), digits = 1),
    across(
      c(starts_with("fecha_kings_"), starts_with("fecha_mitos_")),
      ~ round((.x - fecha_inicio_clinica) / dmonths(1), digits = 1),
      .names = "tiempo_hasta_{str_remove(.col, 'fecha_')}"
    ),
  ), by = "nhc") |>
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

als_samples_info |>
  dplyr::arrange(factor(grupo, levels = c("ELA", "CONTROL")), sample_id) |>
  transmute(
    SAMPLE = sample_id,
    NHC = nhc,
    SEX = sexo,
    PHENO = grupo |> case_match("CONTROL" ~ "CONTROL", "ELA" ~ "ALS"),
    ALS_ONSET = edad_inicio,
    ALS_PHENO = fenotipo,
    ALS_DD = retraso_diagnostico_meses,
    ALS_DFS = delta_fs,
    ALS_PC = categoria_progresion,
    ALS_FVC = fvc_basal,
    ALS_KSS1 = tiempo_hasta_kings_1,
    ALS_KSS2 = tiempo_hasta_kings_2,
    ALS_KSS3 = tiempo_hasta_kings_3,
    ALS_KSS4 = tiempo_hasta_kings_4,
    ALS_MITOS1 = tiempo_hasta_mitos_1,
    ALS_MITOS2 = tiempo_hasta_mitos_2,
    ALS_MITOS3 = tiempo_hasta_mitos_3,
    ALS_MITOS4 = tiempo_hasta_mitos_4,
    ALS_SURVIVAL = supervivencia_meses,
  ) |>
  write_tsv(stdout())
