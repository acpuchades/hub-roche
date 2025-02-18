library(dplyr)
library(tidyr)
library(lubridate)
library(DBI)
library(RSQLite)

source("src/utils.r")

as_genetic_status <- function(x) {
  factor(x, levels = c("Normal", "Alterado"))
}

as_kings_stage <- function(x) {
  ordered(x, levels = 0:5)
}

as_mitos <- function(x) {
  ordered(x, levels = 0:5)
}

as_progression_category <- function(x) {
  ordered(x, levels = c("SP", "NP", "FP"))
}

as_clinical_phenotype <- function(x) {
  dplyr::case_match(
    x,
    c("ELA Espinal", "ELA Bulbar", "ELA Respiratoria") ~ "ALS",
    "Atrofia Muscular Progresiva (AMP)" ~ "PMA",
    "Esclerosis Lateral Primaria (ELP)" ~ "PLS",
    "Parálisis bulbar progresiva" ~ "PBP",
    "Flail arm" ~ "Flail-Arm",
    "Flail leg" ~ "Flail-Leg",
    .default = "Other"
  )
}

as_site_of_onset <- function(x) {
  dplyr::case_match(
    x,
    "ELA Respiratoria" ~ "Respiratory",
    "Hemipléjica (Mills)" ~ "Spinal",
    "Pseudopolineurítica" ~ "Spinal",
    "Atrofia Muscular Progresiva (AMP)" ~ "Spinal",
    c("ELA Espinal", "Flail arm", "Flail leg") ~ "Spinal",
    c("ELA Bulbar", "Parálisis bulbar progresiva") ~ "Bulbar",
  )
}

ufela_db <- dbConnect(SQLite(), "data/ufela/formulario_2023-11-15.sqlite")

ufela_pacientes <- dbReadTable(ufela_db, "pacientes") |>
  rename_with(normalize_names) |>
  mutate(
    across(everything(), \(x) na_if(x, "NS/NC")),
    across(starts_with("fecha_"), lubridate::dmy),
    nhc = as.integer(nhc),
  )

ufela_clinico <- dbReadTable(ufela_db, "datos_clinicos") |>
  rename_with(normalize_names) |>
  mutate(
    across(everything(), \(x) na_if(x, "NS/NC")),
    across(starts_with("fecha_"), lubridate::dmy),
    estado_c9 = as_genetic_status(resultado_estudio_c9),
    estado_sod1 = as_genetic_status(resultado_estudio_sod1),
    retraso_diagnostico_meses = round(
      (fecha_diagnostico_ela - fecha_inicio_clinica) / dmonths(1),
      digits = 1
    )
  )

ufela_nutri <- dbReadTable(ufela_db, "datos_antro") |>
  rename_with(normalize_names) |>
  rename(
    fecha_visita = fecha_visita_datos_antro,
    altura_cm = estatura,
  ) |>
  select(-imc_actual) |>
  group_by(pid) |>
  arrange(fecha_visita) |>
  fill(peso_premorbido, .direction = "down") |>
  ungroup() |>
  mutate(
    across(everything(), \(x) na_if(x, "NS/NC")),
    across(starts_with("fecha_"), lubridate::dmy),
    across(c(indicacion_peg, portador_peg), ~ case_match(.x, "Sí" ~ TRUE, "No" ~ FALSE)),
    across(c(peso, peso_premorbido, altura_cm), as.numeric),
    imc = round(peso / (altura_cm / 100)^2, digits = 1),
    perdida_ponderal = (peso - peso_premorbido) / peso_premorbido,
  ) |>
  select(
    pid, fecha_visita, peso, peso_premorbido, altura_cm, imc,
    perdida_ponderal, indicacion_peg, portador_peg
  )

ufela_respi <- dbReadTable(ufela_db, "fun_res") |>
  rename_with(normalize_names) |>
  rename(
    fecha_visita = fecha_visita_fun_res,
    fvc_sentado_abs = fvc_sentado_absoluto,
    fvc_estirado_abs = fvc_estirado_absoluto,
  ) |>
  mutate(
    across(everything(), \(x) na_if(x, "NS/NC")),
    across(starts_with("fecha_"), lubridate::dmy),
    across(c(starts_with("fvc_"), pcf, pim, pem, pns), as.numeric),
    fvc = coalesce(fvc_sentado, fvc_estirado),
    fvc_abs = coalesce(fvc_sentado_abs, fvc_estirado_abs),
  ) |>
  select(pid, fecha_visita, fvc, fvc_abs, pcf, pim, pem, pns)

ufela_alsfrs <- dbReadTable(ufela_db, "esc_val_ela") |>
  rename_with(normalize_names) |>
  mutate(
    across(everything(), \(x) na_if(x, "NS/NC")),
    across(starts_with("fecha_"), lubridate::dmy)
  ) |>
  left_join(
    ufela_clinico |> select(pid, fecha_diagnostico_ela, fecha_inicio_clinica),
    by = "pid"
  ) |>
  rename(fecha_visita = fecha_visita_esc_val_ela) |>
  mutate(
    fecha_visita = ymd(fecha_visita),
    across(lenguaje:total, \(x) if_else(x != "NS/NC", as.integer(x), NA_integer_)),
    total_bulbar = lenguaje + salivacion + deglucion,
    total_motor_grosero = cama + caminar + subir_escaleras,
    total_respiratorio = disnea + ortopnea + insuficiencia_respiratoria,
    total_motor_fino = total - total_bulbar - total_motor_grosero - total_respiratorio,
    total = if_else(total |> between(0, 48), total, NA_integer_),
    meses_desde_inicio = (fecha_visita - fecha_inicio_clinica) / dmonths(1),
    meses_desde_diagnostico = (fecha_visita - fecha_diagnostico_ela) / dmonths(1),
    delta_fs = round((48 - total) / meses_desde_inicio, digits = 1),
    mitos = as_mitos({
      movement <- (caminar <= 1) | (vestido <= 1)
      swallowing <- deglucion <= 1
      communication <- (lenguaje <= 1) & (escritura <= 1)
      breathing <- (disnea <= 1) | (insuficiencia_respiratoria <= 2)
      movement + swallowing + communication + breathing
    })
  )

ufela_kings <- ufela_alsfrs |>
  left_join(
    ufela_nutri |>
      filter(indicacion_peg) |>
      slice_min(fecha_visita, by = pid, n = 1) |>
      select(pid, fecha_indicacion_peg = "fecha_visita"),
    by = "pid"
  ) |>
  mutate(
    kings_c = as_kings_stage(case_when(
      !is.na(fecha_indicacion_peg) & fecha_visita >= fecha_indicacion_peg ~ 4,
      cortar_con_peg > 0 & cortar_sin_peg > 0 ~ NA,
      cortar_con_peg == 0 & cortar_sin_peg == 0 ~ NA,
      cortar_con_peg > 0 & cortar_sin_peg == 0 ~ 4,
      (disnea == 0) | (insuficiencia_respiratoria < 4) ~ 4,
      TRUE ~ {
        bulbar_region <- total_bulbar < 12
        upper_region <- (escritura < 4) | (cortar_sin_peg < 4)
        lower_region <- caminar < 4
        bulbar_region + upper_region + lower_region
      }
    )),
    kings = coalesce(as_kings_stage(as.integer(kings)), kings_c),
  ) |>
  select(pid, fecha_kings = "fecha_visita", kings)

dbDisconnect(ufela_db)

ufela_mitos_dates <- ufela_pacientes |>
  select(pid) |>
  cross_join(tibble(mitos = as_mitos(1:4))) |>
  bind_rows(ufela_alsfrs |> select(pid, fecha_mitos = "fecha_visita", mitos)) |>
  drop_na(mitos) |>
  filter(mitos > 0) |>
  slice_min(fecha_mitos, by = c(pid, mitos), n = 1, with_ties = FALSE) |>
  group_by(pid) |>
  arrange(mitos, .by_group = TRUE) |>
  fill(fecha_mitos, .direction = "up") |>
  ungroup() |>
  pivot_wider(names_from = mitos, values_from = fecha_mitos, names_prefix = "fecha_mitos_")

ufela_kings_dates <- ufela_pacientes |>
  select(pid) |>
  cross_join(tibble(kings = as_kings_stage(1:4))) |>
  bind_rows(ufela_kings) |>
  drop_na(kings) |>
  filter(kings > 0) |>
  slice_min(fecha_kings, by = c(pid, kings), n = 1, with_ties = FALSE) |>
  group_by(pid) |>
  arrange(kings, .by_group = TRUE) |>
  fill(fecha_kings, .direction = "up") |>
  ungroup() |>
  pivot_wider(names_from = kings, values_from = fecha_kings, names_prefix = "fecha_kings_")

ufela_deltafs_basal <- ufela_alsfrs |>
  drop_na(delta_fs) |>
  slice_min(fecha_visita, by = pid, n = 1, with_ties = FALSE, na_rm = TRUE) |>
  select(pid, fecha_deltafs = "fecha_visita", delta_fs) |>
  mutate(categoria_progresion = as_progression_category(case_when(
    delta_fs < 0.5 ~ "SP",
    delta_fs |> between(0.5, 1) ~ "NP",
    delta_fs > 1 ~ "FP",
  )))

ufela_biometria_basal <- ufela_nutri |>
  drop_na(peso) |>
  slice_min(fecha_visita, by = pid, n = 1, with_ties = FALSE, na_rm = TRUE) |>
  select(pid, fecha_imc = "fecha_visita", altura_cm, peso_basal = "peso", peso_premorbido, imc_basal = "imc", perdida_ponderal)

ufela_fvc_basal <- ufela_respi |>
  drop_na(fvc) |>
  slice_min(fecha_visita, by = pid, n = 1, with_ties = FALSE, na_rm = TRUE) |>
  select(pid, fecha_fvc = "fecha_visita", fvc_basal = "fvc", fvc_basal_abs = "fvc_abs")

ufela_datos <- ufela_pacientes |>
  left_join(ufela_clinico, by = "pid") |>
  mutate(
    edad_inicio = floor((fecha_inicio_clinica - fecha_nacimiento) / dyears(1)),
    edad_diagnostico = floor((fecha_diagnostico_ela - fecha_nacimiento) / dyears(1)),
  ) |>
  left_join(ufela_fvc_basal, by = "pid") |>
  left_join(ufela_deltafs_basal, by = "pid") |>
  left_join(ufela_biometria_basal, by = "pid") |>
  left_join(ufela_kings_dates, by = "pid") |>
  left_join(ufela_mitos_dates, by = "pid")
