# Rscript R/lalonde_collect.R <results_dir> [<csv>]
# Stacks the per-estimate files of R/lalonde_run.R (or, with a second
# argument, R/lalonde_units_run.R) into data/lalonde_att.csv (or <csv>).
suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) > 1) args[2] else "data/lalonde_att.csv"
source("R/runs.R")
read_run(args[1]) |>
  list_rbind() |>
  write_csv(out)
