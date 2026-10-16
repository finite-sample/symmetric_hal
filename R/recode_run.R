# Rscript R/recode_run.R <out_dir> <recodings> <cores>
# Coding sensitivity of HAL predictions on real data. Each recoding negates
# every covariate independently with probability 1/2 (recoding 1 is the
# original coding). HAL and symmetric HAL are refit on each recoding with the
# same CV folds and test split, so the coding is the only thing that changes;
# each fit saves its test-set predictions.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
source(file.path("R", "basis.R"))
source(file.path("R", "nuisance.R"))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
n_recode <- as.integer(args[2])
cores <- as.integer(args[3])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load_data <- function(name) {
  switch(name,
    boston = {
      data(BostonHousing2, package = "mlbench", envir = environment())
      d <- BostonHousing2
      list(y = d$cmedv,
           x = d |>
             select(crim, zn, indus, chas, nox, rm, age, dis, rad, tax,
                    ptratio, b, lstat) |>
             mutate(chas = as.numeric(as.character(chas))))
    },
    ozone = {
      data(Ozone, package = "mlbench", envir = environment())
      d <- na.omit(Ozone)
      list(y = d$V4, x = select(d, V5:V13))
    }
  )
}

flip <- function(x, signs) sweep(as.matrix(x), 2, signs, `*`)

run_regression <- function(name, signs, method) {
  dat <- load_data(name)
  set.seed(1)
  test <- sample(length(dat$y), round(0.3 * length(dat$y)))
  foldid <- sample(rep_len(1:10, length(dat$y) - length(test)))
  x <- flip(dat$x, signs)
  pred <- fit_nuisance(method, x[-test, ], dat$y[-test], "gaussian", foldid,
                       list(x[test, , drop = FALSE]))[[1]]
  tibble(data = name, method = method, row = test, pred = pred,
         y = dat$y[test])
}

d_of <- c(boston = 13, ozone = 9)
set.seed(2026)
patterns <- map(d_of, \(d) {
  c(list(rep(1, d)),
    map(seq_len(n_recode - 1), \(i) sample(c(-1, 1), d, TRUE)))
})
jobs <- expand_grid(data = names(d_of), recode = seq_len(n_recode),
                    method = c("hal", "shal")) |>
  mutate(path = file.path(out_dir, sprintf("%s_%s_%02d.rds", data, method,
                                           recode))) |>
  filter(!file.exists(path))
cat(nrow(jobs), "fits to run\n")

invisible(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  signs <- patterns[[job$data]][[job$recode]]
  start <- proc.time()[["elapsed"]]
  res <- run_regression(job$data, signs, job$method) |>
    mutate(recode = job$recode, n_flipped = sum(signs < 0),
           seconds = proc.time()[["elapsed"]] - start)
  saveRDS(res, job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
