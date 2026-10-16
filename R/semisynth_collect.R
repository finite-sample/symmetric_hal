# Rscript R/semisynth_collect.R <results_dir>
# Stacks the per-estimate files of R/semisynth_run.R into one CSV in data/.
suppressMessages(library(tidyverse))
list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "\\.rds$",
           full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  write_csv("data/semisynth_att.csv")
