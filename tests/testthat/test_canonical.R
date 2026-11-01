test_that("canonical orientation respects signs, units, and row order", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  x <- cbind(c(-2, -1, 0, 1, 5), c(-2, -1, 0, 1, 2))
  signs <- canonical_signs(x)
  expect_equal(unname(signs[, 1]), c(1, 1))
  expect_setequal(signs[, 2], c(1, -1))
  coding <- c(-1, 1)
  changed <- sweep(sweep(x, 2, c(-12, 0.001), `*`), 2, c(20, 5), `+`)
  oriented <- canonical_signs(changed)
  expect_equal(unname(oriented[, 1]), unname(signs[, 1] * coding[1]))
  expect_setequal(oriented[, 2], signs[, 2])
  expect_equal(canonical_signs(x[5:1, ]), signs)
  balanced <- matrix(-2:2, ncol = 1)
  expect_equal(nrow(canonical_signs(balanced * 0.001 + 100)), 2)
})

test_that("ambiguous orientations average without dropping signal", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  x <- matrix(c(-1, 0, 1), ncol = 1)
  y <- c(1, 0, 1)
  fit <- fit_canonical_hal(x, y, 0.05)
  newx <- matrix(c(-0.25, 0.25), ncol = 1)
  expect_equal(nrow(fit$signs), 2)
  expect_equal(fit$predict(newx)[1], fit$predict(newx)[2])
  expect_gt(diff(range(fit$predict(x))), 0.1)
  expect_equal(fit$predict(newx), fit_canonical_hal(-x, y, 0.05)$predict(-newx))
  x <- matrix(-2:2, ncol = 1)
  y <- c(1, 0, 0, 0, 1)
  newx <- matrix(c(-1.5, -0.5, 0.5, 1.5), ncol = 1)
  original <- fit_canonical_hal(x, y, 0.05)
  shifted <- fit_canonical_hal(x * 0.001 + 100, y, 0.05)
  expect_equal(original$predict(newx), shifted$predict(newx * 0.001 + 100),
               tolerance = 1e-7)
})

test_that("canonical CV orients inside each training fold and is equivariant", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  x <- toy_x(30, 2, 701)^3
  newx <- toy_x(16, 2, 702)
  folds <- rep(1:3, 10)
  signs <- c(-1, 1)
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") x[, 1] * x[, 2] else rep(c(0, 1), 15)
    fit <- cv_canonical_hal(x, y, c(0.01, 0.05), folds, family = family)
    changed <- cv_canonical_hal(sweep(x, 2, signs, `*`), y, c(0.01, 0.05),
                                folds, family = family)
    expect_equal(fit$lambda, changed$lambda)
    expect_equal(fit$cv$loss, changed$cv$loss, tolerance = 1e-7)
    expect_equal(fit$predict(newx),
                 changed$predict(sweep(newx, 2, signs, `*`)), tolerance = 1e-7)
    for (k in 1:3) {
      expect_equal(fit$cv_signs[[as.character(k)]],
                   canonical_signs(x[folds != k, ]))
    }
    permutation <- 30:1
    permuted <- cv_canonical_hal(x[permutation, ], y[permutation],
                                 c(0.01, 0.05), folds[permutation],
                                 family = family)
    expect_equal(fit$predict(newx), permuted$predict(newx), tolerance = 1e-7)
  }
})

test_that("canonical nuisance fits restore Gaussian outcome units", {
  source(file.path(root, "R", "canonical.R"), local = TRUE)
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  x <- toy_x(30, 2, 801)
  y <- 8000 + 2000 * x[, 1]^2
  newx <- toy_x(12, 2, 802)
  folds <- rep(1:3, 10)
  result <- fit_nuisance("canonical", x, y, "gaussian", folds, list(newx))
  changed <- fit_nuisance("canonical", -x, 2 * y + 100, "gaussian", folds,
                          list(-newx))
  expect_equal(changed[[1]], 2 * result[[1]] + 100, tolerance = 1e-7)
  certificates <- as.matrix(attr(result, "diagnostics")$cv_certificates)
  expect_true(all(is.finite(certificates)))
})
