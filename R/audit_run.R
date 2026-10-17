# Rscript R/audit_run.R <out_dir> <cores>
# Which common regression estimators change when the covariates are recoded?
# Boston housing, the 70/30 split and 10 CV folds of R/recode_run.R. Every
# estimator is refit, with practitioners' default settings, under recodings
# that carry the same information: the continuous covariates negated, negated
# about their maximum ("years to retirement" for age), shifted, put in other
# units, and the reference level of the one categorical covariate changed.
# "reseed" refits the original coding with another seed, the scale of
# movement a randomized estimator shows anyway. Each fit saves its test-set
# predictions.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
for (f in c("basis", "nuisance")) source(file.path("R", paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
cores <- as.integer(args[2])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data(BostonHousing2, package = "mlbench")
continuous <- c("crim", "zn", "indus", "nox", "rm", "age", "dis", "tax",
                "ptratio", "b", "lstat")
base <- BostonHousing2 |>
  transmute(across(all_of(continuous)),
            chas = as.numeric(as.character(chas)),
            rad = factor(rad), y = cmedv)
set.seed(1)
test <- sample(nrow(base), round(0.3 * nrow(base)))
foldid <- sample(rep_len(1:10, nrow(base) - length(test)))

# A different change of units for each covariate (as dollars to thousands,
# years to months): a common factor would leave many fits unchanged that a
# mix of units does not.
units <- rep(c(100, 0.01, 12), length.out = length(continuous))
codings <- list(
  original = identity,
  negated = \(d) mutate(d, across(all_of(continuous), \(v) -v)),
  negated_from_max = \(d) {
    mutate(d, across(all_of(continuous), \(v) max(v) - v))
  },
  shifted = \(d) mutate(d, across(all_of(continuous), \(v) v + 10 * sd(v))),
  rescaled = \(d) {
    d[continuous] <- map2(d[continuous], units, `*`)
    d
  },
  reference_level = \(d) mutate(d, rad = fct_relevel(rad, levels(rad)[9])),
  reseed = identity
)

# Dummy coding with treatment contrasts, as model.matrix gives by default.
dummies <- function(d, formula = ~ .) {
  model.matrix(formula, select(d, -y))[, -1, drop = FALSE]
}
glmnet_fit <- function(x, y, newx, alpha) {
  fit <- glmnet::cv.glmnet(x, y, alpha = alpha, foldid = foldid)
  as.vector(predict(fit, newx, s = "lambda.min"))
}

fitters <- list(
  ols = \(tr, te, seed) predict(lm(y ~ ., tr), te),
  ols_interactions = \(tr, te, seed) predict(lm(y ~ .^2, tr), te),
  ridge = \(tr, te, seed) glmnet_fit(dummies(tr), tr$y, dummies(te), 0),
  lasso = \(tr, te, seed) glmnet_fit(dummies(tr), tr$y, dummies(te), 1),
  elastic_net = \(tr, te, seed) {
    glmnet_fit(dummies(tr), tr$y, dummies(te), 0.5)
  },
  lasso_interactions = \(tr, te, seed) {
    glmnet_fit(dummies(tr, ~ .^2), tr$y, dummies(te, ~ .^2), 1)
  },
  gam = \(tr, te, seed) {
    f <- reformulate(c(sprintf("s(%s)", continuous), "chas", "rad"), "y")
    as.vector(predict(mgcv::gam(f, data = tr), te))
  },
  mars = \(tr, te, seed) {
    as.vector(predict(earth::earth(y ~ ., tr, degree = 2), te))
  },
  knn = \(tr, te, seed) {
    x <- dummies(tr)
    centre <- colMeans(x)
    spread <- apply(x, 2, sd)
    FNN::knn.reg(scale(x, centre, spread), scale(dummies(te), centre, spread),
                 tr$y, k = 10)$pred
  },
  svm = \(tr, te, seed) {
    as.vector(predict(e1071::svm(dummies(tr), tr$y), dummies(te)))
  },
  cart = \(tr, te, seed) predict(rpart::rpart(y ~ ., tr), te),
  random_forest = \(tr, te, seed) {
    predict(ranger::ranger(y ~ ., tr, seed = seed, num.threads = 1),
            te)$predictions
  },
  boosting = \(tr, te, seed) {
    set.seed(seed)
    fit <- xgboost::xgb.train(
      params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4,
                    nthread = 1, seed = seed),
      data = xgboost::xgb.DMatrix(dummies(tr), label = tr$y), nrounds = 200
    )
    predict(fit, xgboost::xgb.DMatrix(dummies(te)))
  },
  bart = \(tr, te, seed) {
    fit <- dbarts::bart(select(tr, -y), tr$y, select(te, -y), verbose = FALSE,
                        seed = seed, nthread = 1)
    fit$yhat.test.mean
  }
)
for (m in c("hal", "shal", "hal9001", "shal9001")) {
  fitters[[m]] <- local({
    method <- m
    \(tr, te, seed) {
      set.seed(seed)
      fit_nuisance(method, dummies(tr), tr$y, "gaussian", foldid,
                   list(dummies(te)))[[1]]
    }
  })
}

jobs <- expand_grid(method = names(fitters), coding = names(codings)) |>
  mutate(path = file.path(out_dir, sprintf("%s_%s.rds", method, coding))) |>
  filter(!file.exists(path))
cat(nrow(jobs), "fits to run\n")

invisible(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  d <- codings[[job$coding]](base)
  seed <- if (job$coding == "reseed") 2 else 1
  start <- proc.time()[["elapsed"]]
  pred <- tryCatch(fitters[[job$method]](d[-test, ], d[test, ], seed),
                   error = \(e) {
                     message(job$method, " ", job$coding, ": ",
                             conditionMessage(e))
                     NA_real_
                   })
  tibble(method = job$method, coding = job$coding, row = test,
         pred = as.vector(pred), y = base$y[test],
         seconds = proc.time()[["elapsed"]] - start) |>
    saveRDS(job$path)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
