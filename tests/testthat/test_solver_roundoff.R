test_that("active QP constraints remain exact zeros in the certificate", {
  set.seed(51422)
  correlation <- matrix(c(1, 0.5, 0.5, 1), 2)
  x <- pnorm(matrix(rnorm(600), 300, 2) %*% chol(correlation))
  y <- sin(2 * pi * x[, 1]) + 4 * (x[, 2] - 0.5)^2 +
    rnorm(300, sd = sqrt(0.58747968708209264 / 0.5))
  folds <- sample(rep_len(1:5, 300))
  z <- (y - mean(y)) / sd(y)
  reflected <- sweep(x[folds != 5, ], 2, c(1, -1), `*`)
  fit <- fit_hal(reflected, z[folds != 5], 0.003, solver_tol = 1e-16)
  expect_lt(fit$diagnostics[["global_secondary_kkt"]], 1e-8)
  expect_lt(fit$diagnostics[["prediction_residual"]], 1e-8)
  expect_lt(fit$diagnostics[["budget_change"]], 1e-8)
})


test_that("discrete duplicated covariates retain prediction equivariance", {
  set.seed(53373)
  x <- matrix(sample(0:4, 300, replace = TRUE) / 4, 100, 3)
  x[, 3] <- x[, 2]
  y <- sin(2 * pi * x[, 1]) + 4 * (x[, 2] - 0.5)^2 +
    rnorm(100, sd = sqrt(0.575 / 2))
  folds <- sample(rep_len(1:5, 100))
  z <- (y - mean(y)) / sd(y)
  newx <- as.matrix(expand.grid((0:4) / 4, (0:4) / 4))
  newx <- cbind(newx, newx[, 2])
  lambda <- c(0.003, 0.01, 0.03, 0.1, 0.3)
  original <- cv_hal(x, z, lambda, folds)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), 3)))
  for (i in seq_len(nrow(signs))) {
    fit <- cv_hal(sweep(x, 2, signs[i, ], `*`), z, lambda, folds)
    movement <- max(abs(original$predict(newx) -
                          fit$predict(sweep(newx, 2, signs[i, ], `*`))))
    expect_lt(movement * sd(y), 1e-6 * max(1, sd(y)))
    expect_equal(fit$lambda, original$lambda)
  }
})

test_that("binomial penalty paths reach the requested endpoint accurately", {
  source(file.path(root, "R", "semisynth.R"), local = TRUE)
  source(file.path(root, "R", "att.R"), local = TRUE)
  dgp <- semisynth_dgp(read.csv(file.path(root, "data", "lalonde_dw.csv")))
  sim <- semisynth_draw(dgp, 1)
  folds <- make_folds(sim$a, seed = 1)
  train <- which(folds$outer != 1)
  train <- train[folds$inner[[1]]$all != 1]
  x <- sim$w[train, ]
  y <- sim$a[train]
  anchor <- fit_hal(x, y, 0.002, symmetric = FALSE, family = "binomial")
  expect_identical(anchor$lambda, 0.002)
  expect_lt(max(anchor$diagnostics), 1e-6)
  symmetric <- fit_hal(x, y, 0.002, family = "binomial")
  signs <- c(-1, -1, 1, 1, 1, 1, -1, -1)
  reflected <- fit_hal(sweep(x, 2, signs, `*`), y, 0.002, family = "binomial")
  newx <- sim$w[folds$outer == 1, ]
  movement <- max(abs(symmetric$predict(newx) -
                        reflected$predict(sweep(newx, 2, signs, `*`))))
  expect_lt(movement, 1e-7)
  expect_lt(max(reflected$diagnostics), 1e-6)
})
