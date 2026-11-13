test_that("Gaussian nuisance fits preserve outcome units and reflection", {
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  x <- toy_x(90, 2, 402)
  y <- 2 * (x[, 1] < -0.2 & x[, 2] > 0.2) + sin(x[, 1])
  folds <- rep(1:3, 30)
  newx <- toy_x(20, 2, 403)
  fit <- hal_cv(x, y, TRUE, "gaussian", folds, list(newx))
  dollars <- hal_cv(x, 2000 + 8000 * y, TRUE, "gaussian", folds, list(newx))
  reflected <- hal_cv(-x, 2000 + 8000 * y, TRUE, "gaussian", folds,
                      list(-newx))
  expect_equal(dollars[[1]], 2000 + 8000 * fit[[1]], tolerance = 1e-7)
  expect_equal(reflected[[1]], dollars[[1]], tolerance = 1e-7)
  expect_equal(attr(fit, "diagnostics")$lambda,
               attr(dollars, "diagnostics")$lambda)
  expect_lt(max(attr(dollars, "diagnostics")$certificate), 1e-6)
})

test_that("first-order canonical fits preserve sign and unit codings", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  set.seed(901)
  x <- matrix(sample(-8:8, 80, replace = TRUE), 40, 2)
  newx <- rbind(c(-4, 2), c(3, -6), c(1, 1))
  folds <- rep(1:4, 10)
  change <- function(z) {
    sweep(sweep(z, 2, c(-8, 0.125), `*`), 2, c(16, -32), `+`)
  }
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") x[, 1]^2 + x[, 2] else rep(c(0, 1), 20)
    fit <- fit_nuisance("hal9001_canonical", x, y, family, folds, list(newx))
    changed <- fit_nuisance("hal9001_canonical", change(x), y, family, folds,
                            list(change(newx)))
    expect_equal(fit[[1]], changed[[1]], tolerance = 1e-7)
    expect_equal(attr(fit, "diagnostics")$lambda,
                 attr(changed, "diagnostics")$lambda, tolerance = 1e-7)
  }
})

test_that("first-order canonical fits average ambiguous directions", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  x <- matrix(rep(-2:2, 8), ncol = 1)
  newx <- matrix(c(-1.5, 0.5, 1.5), ncol = 1)
  folds <- rep(1:4, each = 10)
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") x[, 1]^2 else as.numeric(abs(x[, 1]) == 2)
    original <- hal9001_canonical_cv(x, y, family, folds, list(newx))
    changed <- hal9001_canonical_cv(-8 * x + 16, y, family, folds,
                                    list(-8 * newx + 16))
    expect_equal(nrow(attr(original, "diagnostics")$orientations), 2)
    expect_equal(original[[1]], changed[[1]], tolerance = 1e-7)
  }
})
