# Rscript R/recode_collect.R <results_dir>
# Stacks the per-fit files of R/recode_run.R into
# data/recode_predictions.csv, one row per fit x test point.
suppressMessages(library(tidyverse))
dir.create("data", showWarnings = FALSE)
source("R/runs.R")
read_run(commandArgs(trailingOnly = TRUE)[1]) |>
  list_rbind() |>
  write_csv("data/recode_predictions.csv")
