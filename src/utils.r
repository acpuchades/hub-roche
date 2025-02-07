library(stringr)
library(stringi)

normalize_names <- function(x) {
  x |>
    stri_trans_general("Latin-ASCII") |>
    str_to_lower() |>
    str_replace_all("%", "_pct_") |>
    str_replace_all("[^A-Za-z0-9]", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

normalize_sample_id <- function(x) {
  x |>
    str_replace("-?ADN[0-9]", "") |>
    str_replace("-ambBED", "") |>
    str_replace_all("-", "")
}

biobank_donor_id <- function(x) {
  x |>
    str_replace_all("\\([0-9]+\\)$", "") |>
    str_trim()
}

normalize_biobank_id <- function(x) {
  x |>
    str_replace_all(" ", "") |>
    str_replace_all("-", "")
}
