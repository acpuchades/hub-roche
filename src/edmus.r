library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

edmus_normalize_names <- function(data) {
  data |>
    rename_with(\(x) x |>
      str_to_lower() |>
      str_replace("^([^a-z])", "_\\1") |>
      str_replace_all("(?:[^a-z0-9]|\\s)+", "_") |>
      str_replace("_*$", ""))
}

edmus_load <- function(path) {
  read_tsv(path,
    locale = locale(encoding = "UTF-16"),
    guess_max = 9999, na = c("", "?")
  ) |>
    edmus_normalize_names()
}

edmus_personal <- edmus_load("data/ufem/BCN4-Personal-571-240422_173008-DEN.txt") |>
  mutate(
    wait_and_see = as.logical(wait_and_see),
    across(ends_with("_unknown_date"), as.logical),
    across(c(
      first_exam, date_of_birth, date_consent_form, created, ms_onset,
      last_modified, last_info, last_clinical_assessment,
      starts_with("irreversible_") & !ends_with("_unknown_date"),
      ends_with("_date") & !ends_with("_unknown_date") & -unknown_decease_date
    ), dmy)
  ) |>
  drop_na(patient_id)

edmus_diagnosis <- edmus_load("data/ufem/BCN4-Diagnosis-571-240422_173055-DEN.txt") |>
  mutate(
    across(ms_onset, dmy),
    disease_course = factor(case_match(
      disease_course,
      1:2 ~ "RR",
      3 ~ "SP-NR",
      4 ~ "SP-R",
      5 ~ "PP-NR",
      6:7 ~ "PP-R"
    ))
  )

edmus_clinical <- edmus_load("data/ufem/BCN4-Clinical-571-240422_174116-DEN.txt") |>
  mutate(across(date, ~ as.Date(.x, format = "%d/%m/%Y")))

edmus_episodes <- edmus_load("data/ufem/BCN4-Episodes-220811_121723-DEP.txt") |>
  mutate(across(date, ~ as.Date(.x, format = "%d/%m/%Y")))

edmus_arr <- edmus_diagnosis |>
  select(patient_id, ms_onset) |>
  left_join(edmus_episodes, by = "patient_id", multiple = "all", na_matches = "never") |>
  summarize(
    arr_y1 = sum((date - ms_onset) <= dyears(1), na.rm = TRUE),
    arr_y3 = sum((date - ms_onset) <= dyears(3), na.rm = TRUE) / 3,
    .by = "patient_id"
  ) |>
  mutate(across(c(starts_with("arr_")), ~ replace_na(.x, 0)))
