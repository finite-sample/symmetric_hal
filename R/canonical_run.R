# Rscript R/canonical_run.R <pilot_run_dir> <out_dir> <cores> [jobs]
for (file in c("basis", "fit", "canonical", "pilot", "runs")) {
  source(file.path("R", paste0(file, ".R")))
}
for (package in c("glmnet", "quadprog", "Matrix")) loadNamespace(package)
args <- commandArgs(trailingOnly = TRUE)
reference_dir <- args[1]
out_dir <- args[2]
cores <- as.integer(args[3])
settings_path <- file.path(reference_dir, "configuration.rds")
settings <- readRDS(settings_path)
grid <- settings$grid
if (length(args) >= 4) {
  selected <- as.integer(strsplit(args[4], ",")[[1]])
  stopifnot(all(selected %in% grid$job), !anyDuplicated(selected))
  grid <- grid[match(selected, grid$job), ]
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
paths <- file.path(out_dir, sprintf("job_%04d.rds", grid$job))
register_run(out_dir, paths, configuration = list(
  grid = grid, lambda = settings$lambda,
  reference = tools::md5sum(c(settings_path,
                              file.path(reference_dir, "manifest.rds")))
))
todo <- which(!file.exists(paths))
cat(length(todo), "of", nrow(grid), "canonical fits to run\n")
invisible(parallel::mclapply(todo, function(i) {
  cell <- grid[i, ]
  save_job(paths[i], {
    started <- proc.time()[["elapsed"]]
    input <- settings$inputs[[paste(cell$d, cell$distribution)]]
    variance <- input$variance[[cell$truth]]
    draw <- pilot_draw(cell, variance)
    x <- draw$x
    y <- draw$y
    centre <- mean(y)
    spread <- sd(y)
    z <- (y - centre) / spread
    fixed <- if (cell$knot_rule == "fixed") {
      rep(list(c(0.25, 0.5, 0.75)), cell$d)
    }
    fit <- cv_canonical_hal(x, z, settings$lambda, draw$folds, knots = fixed)
    reflected_knots <- if (!is.null(fixed)) lapply(fixed, `-`)
    reflected <- cv_canonical_hal(-x, z, settings$lambda, draw$folds,
                                  knots = reflected_knots)
    test_prediction <- centre + spread * fit$predict(input$test_x)
    movement <- max(abs(test_prediction - centre -
                          spread * reflected$predict(-input$test_x)))
    if (movement > 1e-6 * max(1, spread)) {
      stop("canonical reflection check failed: ", movement)
    }
    grid_knots <- if (is.null(fixed)) make_knots(x, 3) else fixed
    integration <- pilot_integration(input$test_x, cell$truth,
                                     cell$distribution, grid_knots)
    prediction <- centre + spread * fit$predict(integration$x)
    risk <- sum(integration$weights *
                  (prediction^2 - 2 * prediction * integration$mean +
                     integration$second))
    test_risk <- mean((test_prediction -
                         pilot_truths[[cell$truth]](input$test_x))^2)
    stopifnot(is.finite(risk), risk >= -1e-10, is.finite(test_risk))
    list(cell = cell, foldid = draw$folds,
         results = data.frame(
           method = "canonical", risk = risk, test_risk = test_risk,
           lambda = fit$lambda, seconds = sum(fit$seconds),
           signal_variance = variance, exact_integration = integration$exact,
           selected_anchor = NA_integer_, max_symmetric_movement = NA_real_,
           outcome_sd = spread
         ),
         orientations = nrow(fit$signs), cv_signs = fit$cv_signs, cv = fit$cv,
         diagnostics = fit$diagnostics, cv_diagnostics = fit$cv_diagnostics,
         reflection_movement = movement,
         elapsed_seconds = proc.time()[["elapsed"]] - started)
  })
  cat(cell$job, cell$truth, cell$distribution, "finished\n")
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
