# Rscript R/causal_collect.R <nuisance_dir>
source("R/runs.R")
args <- commandArgs(trailingOnly = TRUE)
runs <- read_run(args[1])
manifest <- readRDS(file.path(args[1], "manifest.rds"))
packages <- manifest$specification$runtime$packages
records <- lapply(runs, function(run) {
  p <- run$predictions
  stopifnot(all(is.finite(as.matrix(p))), !anyDuplicated(p$row))
  data.frame(method = run$method, coding = run$coding, p)
})
reference <- records[[1]][c("row", "fold", "a", "y")]
stopifnot(all(vapply(records, function(p) {
  identical(p[c("row", "fold", "a", "y")], reference)
}, logical(1))))
readr::write_csv(do.call(rbind, records), "data/causal_nuisance.csv")
summary <- do.call(rbind, lapply(runs, function(run) {
  data.frame(method = run$method, coding = run$coding, seconds = run$seconds,
             hal9001_version = packages["hal9001", "Version"],
             glmnet_version = packages["glmnet", "Version"],
             warnings = length(run$warnings),
             warning_messages = paste(unique(run$warnings), collapse = " | "))
}))
readr::write_csv(summary, "data/causal_fits.csv")
writeLines(trimws(c(capture.output(manifest$specification$runtime),
                    capture.output(manifest$session)), which = "right"),
           "data/causal_session.txt")
