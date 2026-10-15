# Rscript R/recode_collect.R <results_dir>
# Stacks the per-fit files of R/recode_run.R into data/recode_predictions.csv
# (test-set predictions) and data/recode_ate.csv (lalonde estimates).
suppressMessages(library(tidyverse))
results <- list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "\\.rds$",
                      full.names = TRUE) |>
  map(readRDS)
dir.create("data", showWarnings = FALSE)
keep(results, \(r) "ate" %in% names(r)) |>
  list_rbind() |>
  write_csv("data/recode_ate.csv")
discard(results, \(r) "ate" %in% names(r)) |>
  list_rbind() |>
  write_csv("data/recode_predictions.csv")
