test_that("first-order reflection changes exact lasso fits on a common span", {
  x <- c(0, 0.5, 1)
  y <- c(0, 0, 0.5)
  left <- cbind(x, pmax(x - 0.5, 0))
  right <- cbind(1 - x, pmax(0.5 - x, 0))
  lambda <- 1 / 216
  beta_left <- c(1 / 18, 5 / 6)
  beta_right <- c(-13 / 18, 1 / 2)
  fit_left <- as.vector(left %*% beta_left)
  fit_right <- as.vector(4 / 9 + right %*% beta_right)
  expect_equal(fit_left, c(0, 1 / 36, 17 / 36))
  expect_equal(fit_right, c(-1 / 36, 1 / 12, 4 / 9))
  expect_equal(mean(fit_left - y), 0, tolerance = 1e-15)
  expect_equal(mean(fit_right - y), 0, tolerance = 1e-15)
  expect_equal(as.vector(crossprod(left, fit_left - y)) / length(y),
               -lambda * sign(beta_left))
  expect_equal(as.vector(crossprod(right, fit_right - y)) / length(y),
               -lambda * sign(beta_right))
  expect_equal(qr(cbind(1, left))$rank, 3)
  expect_equal(qr(cbind(1, right))$rank, 3)
  expect_equal(max(abs(fit_left - fit_right)), 1 / 18)

  grid <- seq(0, 1, length.out = 101)
  hinge <- pmax(grid - 0.5, 0)
  reflected_representation <- 0.5 - (1 - grid) + pmax(0.5 - grid, 0)
  expect_equal(reflected_representation, hinge)
})

test_that("hal9001 first-order atoms match the analytic reflection example", {
  x <- matrix(c(0, 0.5, 1), ncol = 1)
  for (sign in c(1, -1)) {
    coded <- sign * x
    basis <- hal9001::enumerate_basis(coded, max_degree = 1,
                                      smoothness_orders = 1, num_knots = 3)
    cutoffs <- vapply(basis, function(atom) atom$cutoffs, numeric(1))
    expected_cutoffs <- sort(as.vector(coded))
    expect_equal(sort(cutoffs), expected_cutoffs)
    actual <- as.matrix(hal9001::make_design_matrix(coded, basis))
    expected <- vapply(cutoffs, function(knot) {
      pmax(as.vector(coded) - knot, 0)
    }, numeric(nrow(coded)))
    expect_equal(unname(actual), unname(expected))
  }
})
