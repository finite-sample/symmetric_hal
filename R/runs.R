# A run directory belongs to one complete job grid and one set of source inputs.
register_run <- function(out_dir, paths, inputs = c(
  sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
  "DESCRIPTION", "data/lalonde_dw.csv"
), runtime = list(
  R = R.version.string,
  packages = utils::installed.packages()[, c("Package", "Version")]
), configuration = NULL) {
  jobs <- basename(paths)
  if (!length(jobs) || anyDuplicated(jobs)) stop("invalid run job grid")
  if (any(!file.exists(inputs))) stop("run source inputs are missing")
  specification <- list(jobs = jobs, inputs = tools::md5sum(inputs),
                        runtime = runtime, configuration = configuration)
  path <- file.path(out_dir, "manifest.rds")
  if (file.exists(path)) {
    if (!identical(readRDS(path)$specification, specification)) {
      stop("run inputs changed; use a fresh output directory")
    }
  } else {
    if (length(list.files(out_dir, pattern = "\\.rds$"))) {
      stop("checkpoint files have no manifest; use a fresh output directory")
    }
    sources <- setNames(lapply(inputs, function(input) {
      readBin(input, "raw", n = file.info(input)$size)
    }), inputs)
    saveRDS(list(specification = specification, sources = sources,
                 session = utils::sessionInfo()), path)
  }
  invisible(NULL)
}

read_run <- function(out_dir) {
  path <- file.path(out_dir, "manifest.rds")
  if (!file.exists(path)) stop("run manifest is missing")
  manifest <- readRDS(path)$specification
  paths <- file.path(out_dir, manifest$jobs)
  missing <- manifest$jobs[!file.exists(paths)]
  if (length(missing)) {
    stop("incomplete run: ", length(missing), " missing jobs")
  }
  results <- lapply(paths, readRDS)
  failed <- vapply(results, function(r) {
    inherits(r, "try-error") ||
      (!is.null(r[["error"]]) && any(!is.na(r[["error"]]))) ||
      any(vapply(c("pred", "att", "se"), function(column) {
        value <- r[[column]]
        !is.null(value) && any(!is.finite(value))
      }, logical(1)))
  }, logical(1))
  if (any(failed)) {
    stop("failed jobs: ", paste(manifest$jobs[failed], collapse = ", "))
  }
  results
}

save_job <- function(path, expr) {
  result <- tryCatch(force(expr), error = function(e) {
    list(error = conditionMessage(e))
  })
  saveRDS(result, path)
  invisible(NULL)
}
