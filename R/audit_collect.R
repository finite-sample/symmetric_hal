# Rscript R/audit_collect.R <results_dir>
# Stacks the per-fit files of R/audit_run.R into data/audit_predictions.csv,
# one row per estimator x coding x test point.
suppressMessages(library(tidyverse))
dir.create("data", showWarnings = FALSE)
list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "\\.rds$",
           full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  select(-seconds) |>
  write_csv("data/audit_predictions.csv")
