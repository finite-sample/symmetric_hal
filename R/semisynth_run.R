# Rscript R/semisynth_run.R <out_dir> <reps> <cores> [methods] [n]
files <- c("basis", "fit", "canonical", "nuisance", "att", "causal",
           "semisynth", "runs")
for (file in files) {
  source(file.path("R", paste0(file, ".R")))
}
for (package in c("glmnet", "quadprog", "Matrix", "hal9001")) {
  loadNamespace(package)
}
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
reps <- as.integer(args[2])
cores <- as.integer(args[3])
methods <- if (length(args) >= 4) {
  strsplit(args[4], ",")[[1]]
} else {
  c("oracle", "hal", "canonical", "shal", "hal9001", "hal9001_canonical")
}
dgp <- semisynth_dgp(read.csv("data/lalonde_dw.csv"))
n <- if (length(args) >= 5) as.integer(args[5]) else nrow(dgp$population)
codings <- list(original = rep(1, 8),
                continuous_negated = c(-1, -1, 1, 1, 1, 1, -1, -1))
jobs <- expand.grid(replicate = seq_len(reps), method = methods,
                    coding = names(codings), stringsAsFactors = FALSE)
paths <- file.path(out_dir, sprintf("%03d_%s_%s.rds", jobs$replicate,
                                    jobs$method, jobs$coding))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
register_run(out_dir, paths, configuration = list(
  jobs = jobs, n = n, codings = codings, dgp = dgp
))
todo <- which(!file.exists(paths))
cat(length(todo), "of", length(paths), "jobs to run\n")
invisible(parallel::mclapply(todo, function(i) {
  job <- jobs[i, ]
  save_job(paths[i], {
    sim <- semisynth_draw(dgp, job$replicate, n)
    folds <- make_folds(sim$a, seed = job$replicate)
    warnings <- character()
    start <- proc.time()[["elapsed"]]
    predictions <- if (job$method == "oracle") {
      data.frame(row = seq_len(n), fold = folds$outer, a = sim$a, y = sim$y,
                 q0 = sim$truth$true_q0, q1 = sim$truth$true_q1,
                 g = sim$truth$true_g)
    } else {
      withCallingHandlers(
        cross_fit_nuisance(sweep(sim$w, 2, codings[[job$coding]], `*`),
                           sim$a, sim$y, job$method, folds$outer, folds$inner,
                           include_treated = TRUE),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }
    nuisance_seconds <- proc.time()[["elapsed"]] - start
    list(job = job, predictions = predictions, truth = sim$truth,
         effect = sim$effect, bound_se = sim$bound_se,
         nuisance_loss = semisynth_nuisance_loss(predictions, sim$truth),
         folds = folds, warnings = warnings,
         nuisance_seconds = nuisance_seconds,
         seconds = proc.time()[["elapsed"]] - start)
  })
  cat(job$replicate, job$method, job$coding, "finished\n")
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
