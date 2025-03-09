library(cli)
library(nls.multstart)

source("src/ufela.r")

d50_model <- function(t, D50, Dx) {
  48 / (1 + exp((t - D50) / Dx))
}

if (file.exists("output/phenotype-info/ufela-d50.rds")) {
  cli_inform("Loading cached D50 data from {.file output/phenotype-info/ufela-d50.rds}")
  ufela_d50 <- readRDS("output/phenotype-info/ufela-d50.rds")
} else {
  d50_input <- ufela_alsfrs |>
    transmute(pid, t = (fecha_visita - fecha_inicio_clinica) / dmonths(1), total) |>
    drop_na() |>
    filter(total |> between(1, 48)) |>
    filter(n() >= 3, .by = pid)

  ufela_d50 <- NULL
  pbar <- cli_progress_bar("Estimating patients D50", total = n_distinct(d50_input$pid))
  d50_input |>
    group_by(pid) |>
    group_map(~ {
      cli_progress_update(id = pbar)

      model.fit <- nls_multstart(
        total ~ d50_model(t, D50, Dx),
        data = .x,
        start_lower = c(D50 = 3, Dx = 1),
        start_upper = c(D50 = 10 * 12, Dx = 10),
        iter = 500, supp_errors = "Y",
      )

      if (is.null(model.fit)) {
        return()
      }

      ufela_d50 <<- bind_rows(ufela_d50, tibble(
        pid = .y$pid[1],
        D50 = coef(model.fit)[1],
        Dx = coef(model.fit)[2]
      ))
    })
  cli_progress_done(id = pbar)

  cli_inform("Saving D50 data to {.file output/phenotype-info/ufela-d50.rds}")
  saveRDS(ufela_d50, "output/phenotype-info/ufela-d50.rds")
}
