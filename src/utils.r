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