# Rscript R/causal_run.R <out_dir> <cores> [methods]
suppressPackageStartupMessages(library(tidyverse))
for (file in c("basis", "fit", "canonical", "nuisance", "att", "runs")) {
  source(file.path("R", paste0(file, ".R")))
}
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
cores <- as.integer(args[2])
methods <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else "hal9001"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data <- read.csv("data/lalonde_dw.csv")
data <- data[data$group %in% c("treated", "psid"), ]
columns <- c("age", "educ", "black", "hisp", "married", "nodegree",
             "re74", "re75")
x <- as.matrix(data[columns])
codings <- list(original = rep(1, 8),
                continuous_negated = c(-1, -1, 1, 1, 1, 1, -1, -1),
                units = c(12, 1, 1, 1, 1, 1, 0.001, 0.001))
folds <- make_folds(data$treat, seed = 1)
jobs <- expand.grid(method = methods, coding = names(codings),
                    stringsAsFactors = FALSE)
paths <- file.path(out_dir, paste0(jobs$method, "_", jobs$coding, ".rds"))
register_run(out_dir, paths, configuration = list(
  jobs = jobs, codings = codings, folds = folds
))
todo <- which(!file.exists(paths))
cat(length(todo), "nuisance fits to run\n")
invisible(parallel::mclapply(todo, function(i) {
  job <- jobs[i, ]
  warnings <- character()
  save_job(paths[i], {
    start <- proc.time()[["elapsed"]]
    predictions <- withCallingHandlers(
      cross_fit_nuisance(sweep(x, 2, codings[[job$coding]], `*`),
                         data$treat, data$re78, job$method,
                         folds$outer, folds$inner, include_treated = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(method = job$method, coding = job$coding, predictions = predictions,
         folds = folds, seconds = proc.time()[["elapsed"]] - start,
         warnings = warnings, diagnostics = attr(predictions, "diagnostics"))
  })
  cat(job$method, job$coding, "finished\n")
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
