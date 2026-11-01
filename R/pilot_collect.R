# Rscript R/pilot_collect.R <completed_run_dir> [completed_canonical_dir]
suppressPackageStartupMessages(library(dplyr))
source("R/runs.R")
args <- commandArgs(trailingOnly = TRUE)
runs <- read_run(args[1])
settings <- readRDS(file.path(args[1], "configuration.rds"))
stopifnot(length(runs) == nrow(settings$grid))
records <- do.call(rbind, lapply(runs, function(run) {
  stopifnot(nrow(run$results) == 7, !anyDuplicated(run$results$method),
            all(is.finite(run$results$risk)),
            all(is.finite(run$results$test_risk)))
  data.frame(run$cell, run$results, elapsed_seconds = run$elapsed_seconds,
             row.names = NULL)
}))
stopifnot(setequal(records$job, settings$grid$job))
diagnostics <- do.call(rbind, lapply(runs, function(run) {
  data.frame(run$cell, method = "symmetric",
             max_reflection_movement = max(run$reflection_movement),
             orientations = NA_integer_)
}))
if (length(args) >= 2) {
  canonical <- read_run(args[2])
  stopifnot(length(canonical) == length(runs))
  extra <- lapply(seq_along(runs), function(i) {
    other <- canonical[[i]]
    stopifnot(identical(other$cell, runs[[i]]$cell),
              identical(other$foldid, runs[[i]]$foldid),
              identical(other$results$outcome_sd,
                        runs[[i]]$results$outcome_sd[1]),
              identical(other$results$method, "canonical"))
    data.frame(other$cell, other$results,
               elapsed_seconds = other$elapsed_seconds, row.names = NULL)
  })
  records <- rbind(records, do.call(rbind, extra))
  extra_diagnostics <- lapply(canonical, function(run) {
    data.frame(run$cell, method = "canonical",
               max_reflection_movement = run$reflection_movement,
               orientations = run$orientations)
  })
  diagnostics <- rbind(diagnostics, do.call(rbind, extra_diagnostics))
}
readr::write_csv(records, "data/pilot_runs.csv")
readr::write_csv(diagnostics, "data/pilot_diagnostics.csv")
cells <- c("truth", "d", "n", "snr", "distribution", "knot_rule")
summary <- records |>
  group_by(across(all_of(c(cells, "method")))) |>
  summarise(replicates = n(), risk_mcse = sd(risk) / sqrt(n()),
            risk = mean(risk), test_risk = mean(test_risk),
            mean_seconds = mean(seconds),
            max_symmetric_movement = max(max_symmetric_movement),
            .groups = "drop")
readr::write_csv(summary, "data/pilot_summary.csv")
symmetric <- records[records$method == "symmetric", ]
comparators <- c("anchored", "selected", "averaged", "cart", "rf")
if (length(args) >= 2) comparators <- c(comparators, "canonical")
comparisons <- lapply(comparators, function(method) {
  other <- records[records$method == method, ]
  other <- other[match(symmetric$job, other$job), ]
  stopifnot(identical(other$job, symmetric$job),
            identical(other$replicate, symmetric$replicate))
  metric <- if (method %in% c("cart", "rf")) "test_risk" else "risk"
  data.frame(symmetric[c(cells, "replicate", "job")], comparator = method,
             symmetric_risk = symmetric[[metric]],
             comparator_risk = other[[metric]],
             exact_integration = symmetric$exact_integration & metric == "risk")
})
paired <- do.call(rbind, comparisons)
paired$paired_difference <- paired$symmetric_risk - paired$comparator_risk
contrasts <- paired |>
  group_by(across(all_of(c(cells, "comparator", "exact_integration")))) |>
  summarise(
    log_ratio_mcse = sd(symmetric_risk / mean(symmetric_risk) -
                          comparator_risk / mean(comparator_risk)) / sqrt(n()),
    replicates = n(), symmetric_risk = mean(symmetric_risk),
    comparator_risk = mean(comparator_risk),
    difference = mean(paired_difference),
    difference_mcse = sd(paired_difference) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(risk_ratio = symmetric_risk / comparator_risk,
         ratio_lower = exp(log(risk_ratio) -
                             qt(0.975, replicates - 1) * log_ratio_mcse),
         ratio_upper = exp(log(risk_ratio) +
                             qt(0.975, replicates - 1) * log_ratio_mcse),
         lower = difference - qt(0.975, replicates - 1) * difference_mcse,
         upper = difference + qt(0.975, replicates - 1) * difference_mcse)
readr::write_csv(contrasts, "data/pilot_comparisons.csv")
print(contrasts)
