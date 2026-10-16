# Rscript R/semisynth_run.R <out_dir> <reps> <cores>
# Semi-synthetic LaLonde: covariates are resampled from the 185 trainees and
# the PSID-1 comparison group; treatment and earnings are drawn from known
# functions fitted to the real data, so the effect on the treated ($1,800)
# and its efficiency bound are known. Each replicate estimates the effect
# under the original coding and with the continuous covariates negated.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
for (f in c("basis", "nuisance", "att")) source(file.path("R", paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
reps <- seq_len(as.integer(args[2]))
cores <- as.integer(args[3])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

covariates <- c("age", "educ", "black", "hisp", "married", "nodegree",
                "re74", "re75")
base <- read_csv(file.path("data", "lalonde_dw.csv"),
                 show_col_types = FALSE) |>
  filter(group %in% c("treated", "psid")) |>
  mutate(u74 = as.numeric(re74 == 0), u75 = as.numeric(re75 == 0))

# True propensity: Dehejia and Wahba's logit, capped at 0.95. True earnings
# without training: splines in age and earnings plus the zero-earnings
# steps, fit to the comparison group. The effect is a constant $1,800.
g_model <- suppressWarnings(glm(
  treat ~ age + I(age^2) + educ + black + hisp + married + nodegree + re74 +
    re75 + u74 + u75,
  family = binomial, data = base
))
q_model <- lm(re78 ~ splines::ns(age, 3) + educ + black + hisp + married +
                nodegree + splines::ns(re74, 4) + splines::ns(re75, 4) + u74 +
                u75, data = filter(base, treat == 0))
true_g <- function(d) pmin(predict(g_model, d, type = "response"), 0.95)
true_q0 <- function(d) predict(q_model, d)
effect <- 1800
sigma <- sigma(q_model)

codings <- list(original = rep(1, 8),
                continuous_negated = c(-1, -1, 1, 1, 1, 1, -1, -1))
methods <- c("hal9001", "shal9001")

jobs <- expand_grid(rep = reps, coding = names(codings), method = methods) |>
  mutate(path = file.path(out_dir, sprintf("%03d_%s_%s.rds", rep, coding,
                                           method))) |>
  filter(!file.exists(path))
cat(nrow(jobs), "estimates to run\n")

draw <- function(rep) {
  set.seed(rep)
  d <- slice_sample(base, n = nrow(base), replace = TRUE)
  g <- true_g(d)
  a <- rbinom(nrow(d), 1, g)
  list(w = as.matrix(d[covariates]), a = a,
       y = true_q0(d) + effect * a + rnorm(nrow(d), sd = sigma), g = g)
}

invisible(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  sim <- draw(job$rep)
  folds <- make_folds(sim$a, seed = job$rep)
  w <- sweep(sim$w, 2, codings[[job$coding]], `*`)
  start <- proc.time()[["elapsed"]]
  res <- tryCatch(
    cross_fit_att(w, sim$a, sim$y, job$method, folds$outer, folds$inner),
    error = function(e) tibble(error = conditionMessage(e))
  )
  # Efficiency bound for this draw's design: sigma^2 E[g / (1 - g)] / p^2.
  bound_se <- sqrt(sigma^2 * mean(sim$g / (1 - sim$g)) / mean(sim$g)^2 /
                     length(sim$a))
  res |>
    mutate(rep = job$rep, coding = job$coding, method = job$method,
           truth = effect, bound_se = bound_se,
           seconds = proc.time()[["elapsed"]] - start) |>
    saveRDS(job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
