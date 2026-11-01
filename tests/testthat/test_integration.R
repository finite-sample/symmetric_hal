test_that("nuisance wrappers use the public estimator for both families", {
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  x <- toy_x(30, 2, 101)
  foldid <- rep(1:3, 10)
  newx <- toy_x(10, 2, 102)
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") x[, 1] + x[, 2]^2 else rep(c(0, 1), 15)
    lambda <- if (family == "gaussian") {
      sd(y) * c(0.003, 0.01, 0.03, 0.1)
    } else {
      c(0.002, 0.01, 0.05, 0.2)
    }
    for (symmetric in c(FALSE, TRUE)) {
      expected <- cv_hal(x, y, lambda, foldid,
                         symmetric = symmetric, family = family)
      method <- if (symmetric) "shal" else "hal"
      actual <- fit_nuisance(method, x, y, family, foldid, list(newx))[[1]]
      expect_equal(actual, expected$predict(newx))
    }
  }
  expect_error(fit_nuisance("unknown", x, x[, 1], "gaussian", foldid,
                            list(newx)), "unknown nuisance method")
})

test_that("simulation pipeline returns complete finite paired fits", {
  source(file.path(root, "R", "sim_setup.R"), local = TRUE)
  test_x <- toy_x(40, 2, 104)
  fits <- run_replicate("box_mixed", 2, 1, 40, 2, 0, test_x,
                        c("hal", "shal"), num_knots = 1)
  expect_equal(vapply(fits, `[[`, character(1), "method"), c("hal", "shal"))
  expect_true(all(vapply(fits, function(fit) {
    is.finite(fit$risk) && all(is.finite(fit$pred)) && fit$seconds >= 0
  }, logical(1))))
  expect_equal(vapply(fits, `[[`, numeric(1), "rep"), c(1, 1))
})

test_that("ATT estimation rejects nonfinite nuisance predictions", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  fit_nuisance <- function(method, x, y, family, foldid, newx) {
    lapply(newx, function(z) rep(NA_real_, nrow(z)))
  }
  x <- toy_x(30, 2, 105)
  a <- rep(c(0, 1), 15)
  folds <- make_folds(a)
  expect_error(cross_fit_att(x, a, x[, 1], "hal", folds$outer, folds$inner),
               "invalid nuisance predictions")
})
