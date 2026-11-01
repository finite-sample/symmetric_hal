# Rscript R/semisynth_collect.R <completed_run_dir> [additional_completed_dirs]
for (file in c("runs", "att", "causal")) {
  source(file.path("R", paste0(file, ".R")))
}
directories <- commandArgs(trailingOnly = TRUE)
stopifnot(length(directories) > 0)
manifests <- lapply(file.path(directories, "manifest.rds"), readRDS)
configurations <- lapply(manifests, function(manifest) {
  manifest$specification$configuration
})
stopifnot(all(vapply(configurations, function(configuration) {
  all(c("n", "codings", "dgp") %in% names(configuration))
}, logical(1))))
reference_design <- configurations[[1]][c("n", "codings", "dgp")]
stopifnot(all(vapply(configurations, function(configuration) {
  identical(configuration[c("n", "codings", "dgp")], reference_design)
}, logical(1))))
runs <- unlist(lapply(directories, read_run), recursive = FALSE)
jobs <- do.call(rbind, lapply(runs, `[[`, "job"))
stopifnot(!anyDuplicated(jobs[c("replicate", "method", "coding")]))
for (replicate in unique(vapply(runs, function(run) run$job$replicate, 1L))) {
  paired <- Filter(function(run) run$job$replicate == replicate, runs)
  reference <- paired[[1]]
  stopifnot(all(vapply(paired, function(run) {
    all(vapply(c("row", "fold", "a", "y"), function(column) {
      identical(run$predictions[[column]], reference$predictions[[column]])
    }, logical(1))) &&
      identical(run$truth, reference$truth) &&
      identical(run$effect, reference$effect) &&
      identical(run$bound_se, reference$bound_se)
  }, logical(1))))
}
cores <- as.integer(Sys.getenv("SEMISYNTH_CORES", "1"))
stopifnot(length(cores) == 1, !is.na(cores), cores >= 1)
records <- parallel::mclapply(runs, function(run) {
  start <- proc.time()[["elapsed"]]
  estimates <- compare_att(run$predictions, steps = c(1e-5, 1e-6))
  targeting_seconds <- proc.time()[["elapsed"]] - start
  finest <- min(estimates$step, na.rm = TRUE)
  final <- subset(estimates, estimator != "TMLE" | step == finest)
  if (!isTRUE(all(final$valid_targeting))) {
    stop("targeting checks failed for replicate ", run$job$replicate,
         ", method ", run$job$method, ", coding ", run$job$coding)
  }
  clipping <- c("clipped_q0", "clipped_q1", "clipped_g")
  stopifnot(all(clipping %in% names(final)))
  names(final)[match(clipping, names(final))] <- paste0("prepared_", clipping)
  data.frame(run$job, final, truth = run$effect, bound_se = run$bound_se,
             n = nrow(run$predictions), run$nuisance_loss,
             nuisance_warnings = length(run$warnings),
             nuisance_seconds = run$nuisance_seconds,
             targeting_seconds = targeting_seconds,
             seconds = run$nuisance_seconds + targeting_seconds,
             row.names = NULL)
}, mc.cores = cores)
stopifnot(!any(vapply(records, inherits, logical(1), "try-error")))
readr::write_csv(do.call(rbind, records), "data/semisynth_att.csv")
nuisance_sessions <- setNames(lapply(manifests, function(manifest) {
  list(inputs = manifest$specification$inputs, session = manifest$session)
}), basename(directories))
writeLines(trimws(c(
  capture.output(nuisance_sessions),
  capture.output(tools::md5sum(c("R/att.R", "R/causal.R",
                                 "R/semisynth_collect.R"))),
  capture.output(sessionInfo())
), which = "right"), "data/semisynth_session.txt")
