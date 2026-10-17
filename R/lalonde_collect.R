# Rscript R/lalonde_collect.R <results_dir> [<csv>]
# Stacks the per-estimate files of R/lalonde_run.R (or, with a second
# argument, R/lalonde_units_run.R) into data/lalonde_att.csv (or <csv>).
suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) > 1) args[2] else "data/lalonde_att.csv"
list.files(args[1], pattern = "\\.rds$", full.names = TRUE) |>
  map(readRDS) |>
  list_rbind() |>
  write_csv(out)
