# Rscript R/recode_run.R <out_dir> <recodings> <cores>
# Coding sensitivity on real data. Each recoding negates every covariate
# independently with probability 1/2 (recoding 1 is the original coding).
# HAL and symmetric HAL are refit on each recoding with the same CV folds,
# so the coding is the only thing that changes. Regression data sets save
# test-set predictions; lalonde saves a one-step (AIPW) estimate of the
# average treatment effect of training on 1978 earnings and its standard
# error, with HAL fits for the outcome regression and the propensity score.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
source(file.path("R", "basis.R"))

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
    },
    lalonde = {
      data(lalonde, package = "MatchIt", envir = environment())
      d <- lalonde
      list(y = d$re78, a = d$treat,
           x = d |>
             mutate(black = as.numeric(race == "black"),
                    hispan = as.numeric(race == "hispan")) |>
             select(age, educ, black, hispan, married, nodegree, re74, re75))
    }
  )
}

hal_cv <- function(x, y, symmetric, family, foldid, newx) {
  features <- make_dictionary(x, 9, symmetric)
  # The default convergence threshold leaves CV curves loose enough that
  # lambda.min can differ between equivalent codings; tighten it so any
  # remaining difference is the estimator's, not the solver's.
  cv <- glmnet::cv.glmnet(design_matrix(x, features), y, family = family,
                          standardize = FALSE, foldid = foldid,
                          control = list(thresh = 1e-12))
  map(newx, \(nx) {
    as.vector(predict(cv, design_matrix(nx, features), s = "lambda.min",
                      type = "response"))
  })
}

flip <- function(x, signs) sweep(as.matrix(x), 2, signs, `*`)

run_regression <- function(name, signs, method) {
  dat <- load_data(name)
  set.seed(1)
  test <- sample(length(dat$y), round(0.3 * length(dat$y)))
  foldid <- sample(rep_len(1:10, length(dat$y) - length(test)))
  x <- flip(dat$x, signs)
  pred <- hal_cv(x[-test, ], dat$y[-test], method == "shal", "gaussian",
                 foldid, list(x[test, , drop = FALSE]))[[1]]
  tibble(data = name, method = method, row = test, pred = pred,
         y = dat$y[test])
}

run_lalonde <- function(signs, method) {
  dat <- load_data("lalonde")
  n <- length(dat$y)
  set.seed(1)
  foldid <- sample(rep_len(1:10, n))
  w <- flip(dat$x, signs)
  sym <- method == "shal"
  q <- hal_cv(cbind(a = dat$a, w), dat$y, sym, "gaussian", foldid,
              list(cbind(a = dat$a, w), cbind(a = 1, w), cbind(a = 0, w)))
  g <- hal_cv(w, dat$a, sym, "binomial", foldid, list(w))[[1]] |>
    pmin(0.99) |>
    pmax(0.01)
  ic <- q[[2]] - q[[3]] +
    (dat$a / g - (1 - dat$a) / (1 - g)) * (dat$y - q[[1]])
  tibble(data = "lalonde", method = method, ate = mean(ic),
         se = sd(ic) / sqrt(n), min_g = min(g), max_g = max(g))
}

d_of <- c(boston = 13, ozone = 9, lalonde = 8)
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
  res <- if (job$data == "lalonde") run_lalonde(signs, job$method) else
    run_regression(job$data, signs, job$method)
  res <- mutate(res, recode = job$recode,
                n_flipped = sum(signs < 0),
                seconds = proc.time()[["elapsed"]] - start)
  saveRDS(res, job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
