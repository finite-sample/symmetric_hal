# Rscript R/lalonde_run.R <out_dir> <recodings> <cores> [methods]
# LaLonde's test: the effect of NSW training on 1978 earnings for the
# treated, estimated from the 185 Dehejia-Wahba trainees and the PSID-1
# comparison group, against the experimental benchmark. Each estimate is a
# cross-fitted augmented IPW estimate with HAL nuisance fits; it is repeated
# under random sign recodings of the eight covariates (recoding 1 is the
# original coding), with the same outer and inner folds throughout.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
for (f in c("basis", "nuisance", "att")) source(file.path("R", paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
n_recode <- as.integer(args[2])
cores <- as.integer(args[3])
methods <- if (length(args) > 3) strsplit(args[4], ",")[[1]] else
  c("hal", "shal", "hal9001", "shal9001")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

covariates <- c("age", "educ", "black", "hisp", "married", "nodegree",
                "re74", "re75")
dat <- read_csv(file.path("data", "lalonde_dw.csv"), show_col_types = FALSE) |>
  filter(group %in% c("treated", "psid"))
w0 <- as.matrix(dat[covariates])
a <- dat$treat
y <- dat$re78

folds <- make_folds(a, seed = 1)

set.seed(2026)
patterns <- c(list(rep(1, length(covariates))),
              map(seq_len(n_recode - 1), \(i) {
                sample(c(-1, 1), length(covariates), TRUE)
              }))
jobs <- expand_grid(recode = seq_len(n_recode), method = methods) |>
  mutate(path = file.path(out_dir, sprintf("%s_%02d.rds", method, recode))) |>
  filter(!file.exists(path))
cat(nrow(jobs), "estimates to run\n")

invisible(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  signs <- patterns[[job$recode]]
  start <- proc.time()[["elapsed"]]
  res <- cross_fit_att(sweep(w0, 2, signs, `*`), a, y, job$method,
                       folds$outer, folds$inner) |>
    mutate(method = job$method, recode = job$recode,
           flipped = paste(covariates[signs < 0], collapse = " "),
           seconds = proc.time()[["elapsed"]] - start)
  saveRDS(res, job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
