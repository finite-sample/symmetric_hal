# Rscript R/pilot_run.R <out_dir> <feasibility|pilot> <cores> [replicates]
for (file in c("basis", "fit", "anchors", "pilot", "runs")) {
  source(file.path("R", paste0(file, ".R")))
}
for (package in c("glmnet", "quadprog", "Matrix", "rpart", "ranger")) {
  loadNamespace(package)
}
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
phase <- args[2]
cores <- as.integer(args[3])
reps <- if (length(args) >= 4) as.integer(args[4]) else 20L
lambda <- c(0.003, 0.01, 0.03, 0.1, 0.3)
grid <- pilot_grid(reps)
if (phase == "feasibility") {
  grid <- subset(grid, replicate == 1 & n == 100 & snr == 2 &
                   ((d == 2 & distribution == "independent") |
                      (d == 3 & distribution == "discrete" &
                         truth == "conflicting")))
} else if (phase != "pilot") {
  stop("unknown study phase")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
paths <- file.path(out_dir, sprintf("job_%04d.rds", grid$job))
register_run(out_dir, paths, configuration = list(grid = grid, lambda = lambda))
settings_path <- file.path(out_dir, "configuration.rds")
if (!file.exists(settings_path)) {
  designs <- unique(grid[c("d", "distribution")])
  inputs <- lapply(seq_len(nrow(designs)), function(i) {
    design <- designs[i, ]
    seed <- match(design$distribution,
                  c("independent", "correlated", "discrete")) * 10 + design$d
    set.seed(7000 + seed)
    calibration_x <- pilot_x(200000, design$d, design$distribution)
    variance <- vapply(pilot_truths, function(f) {
      var(f(calibration_x))
    }, numeric(1))
    if (design$distribution != "correlated") {
      for (truth in names(pilot_truths)) {
        integration <- pilot_integration(
          matrix(0, 1, design$d), truth, design$distribution,
          rep(list(c(0.25, 0.5, 0.75)), design$d)
        )
        variance[truth] <- sum(integration$weights * integration$second) -
          sum(integration$weights * integration$mean)^2
      }
    }
    set.seed(8000 + seed)
    list(test_x = pilot_x(10000, design$d, design$distribution),
         variance = variance)
  })
  names(inputs) <- paste(designs$d, designs$distribution)
  saveRDS(list(grid = grid, lambda = lambda, inputs = inputs), settings_path)
}
settings <- readRDS(settings_path)
todo <- which(!file.exists(paths))
cat(length(todo), "of", nrow(grid), "replicates to run\n")
invisible(parallel::mclapply(todo, function(i) {
  cell <- grid[i, ]
  input <- settings$inputs[[paste(cell$d, cell$distribution)]]
  save_job(paths[i], {
    pilot_replicate(cell, input$test_x, input$variance[[cell$truth]], lambda)
  })
  cat(cell$job, cell$truth, cell$distribution, "finished\n")
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
