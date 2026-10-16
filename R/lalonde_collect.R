# Rscript R/lalonde_collect.R <results_dir>
# Stacks the per-estimate files of R/lalonde_run.R into data/lalonde_att.csv.
suppressMessages(library(tidyverse))
list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "\\.rds$",
           full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  write_csv("data/lalonde_att.csv")
