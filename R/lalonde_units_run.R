# Rscript R/lalonde_units_run.R <out_dir> <cores>
# The LaLonde estimate of R/lalonde_run.R with the covariates in other units:
# age in months and 1974 and 1975 earnings in thousands of dollars. Same
# sample, folds and estimator; only the units change.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
for (f in c("basis", "nuisance", "att")) source(file.path("R", paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
cores <- as.integer(args[2])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

covariates <- c("age", "educ", "black", "hisp", "married", "nodegree",
                "re74", "re75")
dat <- read_csv(file.path("data", "lalonde_dw.csv"), show_col_types = FALSE) |>
  filter(group %in% c("treated", "psid"))
w0 <- as.matrix(dat[covariates])
a <- dat$treat
y <- dat$re78
folds <- make_folds(a, seed = 1)

units <- list(original = rep(1, length(covariates)),
              months_thousands = c(12, 1, 1, 1, 1, 1, 0.001, 0.001))
jobs <- expand_grid(units = names(units),
                    method = c("hal9001", "shal9001", "shal9001_scaled")) |>
  mutate(path = file.path(out_dir, sprintf("%s_%s.rds", method, units))) |>
  filter(!file.exists(path))
cat(nrow(jobs), "estimates to run\n")

invisible(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  start <- proc.time()[["elapsed"]]
  cross_fit_att(sweep(w0, 2, units[[job$units]], `*`), a, y, job$method,
                folds$outer, folds$inner) |>
    mutate(method = job$method, units = job$units,
           seconds = proc.time()[["elapsed"]] - start) |>
    saveRDS(job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
