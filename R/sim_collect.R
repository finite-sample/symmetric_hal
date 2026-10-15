# Rscript R/sim_collect.R <results_dir> <name>
# Writes data/<name>_runs.csv, one row per replicate x method.
suppressMessages(library(tidyverse))
source(file.path("R", "sim_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
dir.create("data", showWarnings = FALSE)
results <- list.files(args[1], pattern = "\\.rds$", full.names = TRUE) |>
  map(readRDS)

errored <- keep(results, \(r) !is.null(r$error))
if (length(errored) > 0) {
  cat(length(errored), "replicates errored:\n")
  walk(errored, \(r) cat(" ", r$error, "\n"))
}

# Earlier tier-1 replicates also hold methods since dropped from the study.
fits <- discard(results, \(r) !is.null(r$error)) |>
  list_flatten() |>
  map(\(z) {
    as_tibble(z[setdiff(names(z), c("vars", "pred", "parent_cancelled"))])
  }) |>
  list_rbind() |>
  filter(method %in% methods_all)
write_csv(fits, file.path("data", paste0(args[2], "_runs.csv")))
cat(nrow(fits), "runs\n")
