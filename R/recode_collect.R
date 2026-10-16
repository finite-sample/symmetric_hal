# Rscript R/recode_collect.R <results_dir>
# Stacks the per-fit files of R/recode_run.R into
# data/recode_predictions.csv, one row per fit x test point.
suppressMessages(library(tidyverse))
dir.create("data", showWarnings = FALSE)
list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "\\.rds$",
           full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  write_csv("data/recode_predictions.csv")
