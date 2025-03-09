library(tibble)
library(readr)
library(dplyr)
library(readxl)
library(janitor)
library(stringr)
library(tidyr)
library(lubridate)
library(writexl)

source("src/utils.r")
source("src/ufela.r")

discard_samples_cdr <- tibble(codigo_donante = c("ELA959"))
discard_samples_ccn <- tibble(codigo_caso_noraybanks = c("ELA1912", "ELA1923", "ELA1924"))

args <- commandArgs(trailingOnly = TRUE)
input_path <- coalesce(args[1], "output/renamed-variant-samples/samples.normalized.txt")

sample_ids <- read_lines(input_path) |>
  as_tibble_col("sample_id") |>
  filter(str_starts(sample_id, "ELA"))

bb_groups <- read_excel("data/biobanco/Muestras ELA enviadas 20250309.xlsx") |>
  clean_names() |>
  select(-c(x5, x6, x7)) |>
  rename(codigo_caso_noraybanks = codigo_muestra_noray_banks) |>
  mutate(
    nhc = suppressWarnings(as.integer(nhc)),
    across(c(codigo_caso_noraybanks, codigo_donacion_recepcion), normalize_biobank_id),
    codigo_donante = biobank_donor_id(codigo_donacion_recepcion),
    grupo = case_when(
      str_to_upper(tipo_caso) == "CONTROL" ~ "CONTROL",
      str_to_upper(tipo_caso) |> str_detect("ELA") ~ "ELA"
    )
  )

bb_samples <- read_excel("data/biobanco/Muestras ELA.xlsx", skip = 1) |>
  clean_names() |>
  rename(codigo_caso_noraybanks = codigo_caso_noray_banks) |>
  transmute(
    across(nhc, as.integer),
    across(c(codigo_caso_noraybanks, codigo_donacion_recepcion), normalize_biobank_id),
    codigo_donante = biobank_donor_id(codigo_donacion_recepcion),
    fecha_muestra = dmy(coalesce(fecha_muestra, fecha_donacion_recepcion))
  ) |>
  left_join(
    bb_groups |>
      select(codigo_donante, nhc_cdr = nhc, grupo_cdr = grupo) |>
      group_by(codigo_donante) |>
      fill(nhc_cdr, grupo_cdr) |>
      ungroup() |>
      distinct(),
    by = "codigo_donante", na_matches = "never"
  ) |>
  left_join(
    bb_groups |>
      select(codigo_caso_noraybanks, nhc_ccn = nhc, grupo_ccn = grupo) |>
      anti_join(discard_samples_ccn, by = "codigo_caso_noraybanks") |>
      group_by(codigo_caso_noraybanks) |>
      fill(nhc_ccn, grupo_ccn) |>
      ungroup() |>
      distinct(),
    by = "codigo_caso_noraybanks", na_matches = "never"
  ) |>
  left_join(
    ufela_datos |> transmute(nhc, grupo_ufela = "ELA"),
    by = "nhc"
  ) |>
  mutate(
    nhc = coalesce(nhc, nhc_cdr, nhc_ccn),
    grupo = coalesce(grupo_cdr, grupo_ccn, grupo_ufela),
    .keep = "unused"
  )

dropbox_samples <- suppressMessages(read_excel("data/ufela/biobanco-20240201.xlsx")) |>
  rename(codigo_donacion_recepcion = `...2`) |>
  clean_names() |>
  filter(codigo_donacion_recepcion != "Biobanc") |>
  rename(nhc = sap, codigo_caso_noraybanks = caso_noraybanks) |>
  separate_longer_delim(codigo_donacion_recepcion, "/") |>
  mutate(
    nhc = suppressWarnings(as.integer(nhc)),
    across(c(codigo_donacion_recepcion, codigo_caso_noraybanks), normalize_biobank_id),
    codigo_donante = biobank_donor_id(codigo_donacion_recepcion),
    across(starts_with("fecha_"), ~ suppressWarnings(coalesce(
      dmy(.x), as_date(as.numeric(.x), origin = "1899-12-30")
    ))),
    grupo = case_when(
      str_to_lower(dtco) |> str_detect("control") ~ "CONTROL",
      str_to_lower(dtco) |> str_detect("ela|elp|amp|pma|pbp|(flail arm|leg)") ~ "ELA",
      nhc %in% ufela_datos$nhc ~ "ELA",
    )
  )

roche_groups <- bind_rows(
  read_excel("data/biobanco/Muestras-Diagnóstico-20240112.xlsx", col_names = "codigo_donacion_recepcion", sheet = "ELA") |> mutate(grupo = "ELA"),
  read_excel("data/biobanco/Muestras-Diagnóstico-20240112.xlsx", col_names = "codigo_donacion_recepcion", sheet = "Controles") |> mutate(grupo = "Control"),
) |>
  mutate(grupo = grupo |> case_match(
    "Control" ~ "CONTROL",
    "ELA" ~ "ELA"
  ))

samples_info <- sample_ids |>
  left_join(
    bb_samples |>
      select(codigo_caso_noraybanks, nhc_bb_ccn = nhc, grupo_bb_ccn = grupo) |>
      group_by(codigo_caso_noraybanks) |>
      fill(nhc_bb_ccn, grupo_bb_ccn, .direction = "downup") |>
      ungroup() |>
      distinct(),
    by = c(sample_id = "codigo_caso_noraybanks"), na_matches = "never"
  ) |>
  left_join(
    bb_samples |>
      select(codigo_donante, nhc_bb_cdr = nhc, grupo_bb_cdr = grupo) |>
      group_by(codigo_donante) |>
      fill(nhc_bb_cdr, grupo_bb_cdr, .direction = "downup") |>
      ungroup() |>
      distinct(),
    by = c(sample_id = "codigo_donante"), na_matches = "never"
  ) |>
  left_join(
    dropbox_samples |>
      select(codigo_caso_noraybanks, nhc_dropbox_ccn = nhc, grupo_dropbox_ccn = grupo) |>
      anti_join(discard_samples_ccn, by = "codigo_caso_noraybanks") |>
      group_by(codigo_caso_noraybanks) |>
      fill(nhc_dropbox_ccn, grupo_dropbox_ccn, .direction = "downup") |>
      ungroup() |>
      distinct(),
    by = c(sample_id = "codigo_caso_noraybanks"), na_matches = "never"
  ) |>
  left_join(
    dropbox_samples |>
      select(codigo_donante, nhc_dropbox_cdr = nhc, grupo_dropbox_cdr = grupo) |>
      anti_join(discard_samples_cdr, by = "codigo_donante") |>
      group_by(codigo_donante) |>
      fill(nhc_dropbox_cdr, grupo_dropbox_cdr, .direction = "downup") |>
      ungroup() |>
      distinct(),
    by = c(sample_id = "codigo_donante"), na_matches = "never"
  ) |>
  mutate(
    nhc = coalesce(nhc_bb_ccn, nhc_bb_cdr, nhc_dropbox_ccn, nhc_dropbox_cdr),
    grupo = coalesce(grupo_bb_ccn, grupo_bb_cdr, grupo_dropbox_ccn, grupo_dropbox_cdr),
  ) |>
  left_join(
    ufela_datos |>
      transmute(
        nhc,
        sexo = sexo |> case_match(
          "Hombre" ~ "M",
          "Mujer" ~ "F"
        ),
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
      ),
    by = "nhc", na_matches = "never"
  )

samples_info |>
  dplyr::arrange(
    factor(grupo, levels = c("ELA", "CONTROL")),
    sample_id
  ) |>
  transmute(
    SAMPLE = sample_id,
    NHC = nhc,
    SEX = sexo,
    PHENO = grupo |> case_match(
      "CONTROL" ~ "CONTROL",
      "ELA" ~ "ALS"
    ),
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
