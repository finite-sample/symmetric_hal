# Rscript R/sim_collect.R <results_dir> <name>
# Writes data/<name>_runs.csv, one row per replicate x method.
suppressMessages(library(tidyverse))
source(file.path("R", "sim_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
dir.create("data", showWarnings = FALSE)
source("R/runs.R")
results <- read_run(args[1])

fits <- results |>
  list_flatten() |>
  map(\(z) {
    as_tibble(z[setdiff(names(z), c("vars", "pred", "parent_cancelled"))])
  }) |>
  list_rbind() |>
  filter(method %in% methods_all)
write_csv(fits, file.path("data", paste0(args[2], "_runs.csv")))
cat(nrow(fits), "runs\n")
