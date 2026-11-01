reflections <- function(d) {
  as.matrix(expand.grid(rep(list(c(1, -1)), d)))
}

test_that("one-knot anchored fits match the analytic two-group solution", {
  x <- matrix(seq_len(20), ncol = 1)
  y <- rep(c(0, 1), each = 10)
  for (family in c("gaussian", "binomial")) {
    fit <- fit_hal(x, y, 0.1, num_knots = 1, symmetric = FALSE,
                   family = family)
    expect_equal(fit$predict(matrix(c(1, 20), ncol = 1)), c(0.2, 0.8),
                 tolerance = 1e-8)
  }
})

test_that("HAL preserves the lasso optimum and checks both stages", {
  x <- toy_x(40, 2, 41)
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") {
      x[, 1] + x[, 2]^2
    } else {
      as.numeric(x[, 1] + x[, 2] > 0)
    }
    fit <- fit_hal(x, y, 0.04, num_knots = 2, family = family)
    design <- design_matrix(x, fit$features)
    comparison <- glmnet::glmnet(design, y,
      family = family, lambda = 0.04,
      standardize = FALSE,
      control = list(thresh = 1e-16, maxit = 1000000)
    )
    eta <- as.vector(predict(comparison, design, type = "link"))
    beta <- as.vector(coef(comparison))[-1]
    objective <- mean(hal_loss(eta, y, family)) + 0.04 * sum(abs(beta))
    expect_equal(fit$objective, objective, tolerance = 1e-7)
    expect_equal(fit$predict(x, "link"), eta, tolerance = 1e-5)
    expect_lte(sum(fit$beta^2), sum(beta^2) + 1e-6)
    expect_lt(max(fit$diagnostics), 1e-6)
  }
})

test_that("duplicate predictors have equivariant off-diagonal predictions", {
  x <- cbind(1:21, 1:21)
  y <- as.numeric(x[, 1] >= 11)
  newx <- rbind(c(1, 21), c(21, 1), c(5, 17), c(17, 5))
  for (family in c("gaussian", "binomial")) {
    original <- fit_hal(x, y, 0.01, family = family)
    expected <- original$predict(newx)
    expect_equal(expected[1], expected[2], tolerance = 1e-7)
    expect_equal(expected[3], expected[4], tolerance = 1e-7)
    signs <- reflections(2)
    for (i in 2:nrow(signs)) {
      fit <- fit_hal(sweep(x, 2, signs[i, ], `*`), y, 0.01,
        family = family
      )
      actual <- fit$predict(sweep(newx, 2, signs[i, ], `*`))
      expect_equal(actual, expected, tolerance = 1e-6)
      expect_lt(max(fit$diagnostics), 1e-6)
    }
  }
})

test_that("all low-dimensional reflections preserve held-out predictions", {
  for (d in 1:3) {
    continuous <- toy_x(32, d, 70 + d)
    tied <- round(continuous * 2)
    sparse <- tied
    if (d > 1) sparse[, d] <- 1
    scenarios <- list(continuous, tied, sparse)
    for (x in scenarios) {
      y <- x[, 1] + sin(seq_len(nrow(x))) / 3
      newx <- rbind(toy_x(12, d, 80 + d), x[1:5, , drop = FALSE])
      original <- fit_hal(x, y, 0.07, num_knots = 2)
      signs <- reflections(d)
      for (i in 2:nrow(signs)) {
        fit <- fit_hal(sweep(x, 2, signs[i, ], `*`), y, 0.07,
          num_knots = 2
        )
        expect_equal(fit$predict(sweep(newx, 2, signs[i, ], `*`)),
          original$predict(newx),
          tolerance = 1e-6
        )
      }
    }
  }
})

test_that("CV rebuilds training dictionaries and is reflection-equivariant", {
  x <- toy_x(24, 2, 91)
  x[, 2] <- round(x[, 2] * 2)
  foldid <- rep(1:3, 8)
  for (family in c("gaussian", "binomial")) {
    y <- if (family == "gaussian") x[, 1]^2 + x[, 2] else rep(c(0, 1), 12)
    original <- cv_hal(x, y, c(0.08, 0.2), foldid,
      num_knots = 2, family = family
    )
    signs <- reflections(2)
    newx <- toy_x(10, 2, 92)
    for (i in 2:nrow(signs)) {
      fit <- cv_hal(sweep(x, 2, signs[i, ], `*`), y, c(0.2, 0.08),
        foldid,
        num_knots = 2, family = family
      )
      expect_equal(fit$lambda, original$lambda)
      expect_equal(fit$cv$loss, original$cv$loss, tolerance = 1e-6)
      expect_equal(fit$predict(sweep(newx, 2, signs[i, ], `*`)),
        original$predict(newx),
        tolerance = 1e-6
      )
    }
    independent_loss <- numeric(nrow(x))
    for (fold in 1:3) {
      train <- foldid != fold
      fit <- fit_hal(x[train, ], y[train], original$lambda,
        num_knots = 2, family = family
      )
      independent_loss[!train] <- hal_loss(
        fit$predict(x[!train, ], "link"), y[!train], family
      )
    }
    expect_equal(
      mean(independent_loss),
      original$cv$loss[original$cv$lambda == original$lambda]
    )
  }
})

test_that("constant signals and degenerate outcomes have explicit behavior", {
  x <- matrix(rep(1, 24), ncol = 1)
  fit <- cv_hal(x, rep(3, 24), c(0.01, 1, 0.2), rep(1:3, 8))
  expect_equal(fit$lambda, 1)
  expect_equal(fit$predict(matrix(c(-5, 1, 10), ncol = 1)), rep(3, 3),
    tolerance = 1e-7
  )
  expect_error(
    fit_hal(x, rep(0, 24), 0.1, family = "binomial"),
    "both 0 and 1"
  )
  expect_error(
    cv_hal(x, rep(c(0, 1), each = 12), 0.1,
      rep(1:2, each = 12),
      family = "binomial"
    ),
    "both 0 and 1"
  )
  expect_error(fit_hal(x, rep(3, 24), 0), "positive")
  expect_error(fit_hal(x, rep(3, 24), 0.1, num_knots = 0), "integer")
  expect_error(cv_hal(x, rep(3, 24), 0.1, rep(1, 24)), "two explicit")
})

test_that("face selection accounts for primary error on larger responses", {
  x <- cbind(1:21, 1:21)
  y <- 100 * (sin((1:21) / 2) + ((1:21) >= 11))
  newx <- as.matrix(expand.grid(seq(0, 22, by = 2), seq(0, 22, by = 2)))
  original <- fit_hal(x, y, 0.01)
  signs <- reflections(2)
  for (i in 2:nrow(signs)) {
    fit <- fit_hal(sweep(x, 2, signs[i, ], `*`), y, 0.01)
    movement <- max(abs(fit$predict(sweep(newx, 2, signs[i, ], `*`)) -
                          original$predict(newx))) / sd(y)
    expect_lt(movement, 1e-6)
    expect_lt(fit$diagnostics[["global_secondary_kkt"]], 1e-6)
  }
})
