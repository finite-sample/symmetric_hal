and_upper_data <- function(n = 200, seed = 4) {
  x <- toy_x(n = n, d = 3, seed = seed)
  set.seed(seed + 1)
  list(x = x, y = 1.5 * (x[, 1] >= 0.2 & x[, 2] >= 0.2) + rnorm(n, sd = 0.3))
}

# Both dictionaries run through hal9001 unchanged on the design cbind(x, -x),
# so symmetric HAL inherits its families, CV and prediction.
test_that("our lasso matches hal9001::fit_hal at a fixed lambda", {
  dat <- and_upper_data()
  newx <- toy_x(n = 50, d = 3, seed = 9)
  for (sym in c(FALSE, TRUE)) {
    f <- make_dictionary(dat$x, 5, sym)
    ours <- glmnet::glmnet(design_matrix(dat$x, f), dat$y, lambda = 0.01,
                           standardize = FALSE,
                           control = list(thresh = 1e-14))
    theirs <- hal9001::fit_hal(
      cbind(dat$x, -dat$x), dat$y, basis_list = to_hal9001(f, 3),
      smoothness_orders = 0, max_degree = 2, lambda = 0.01,
      reduce_basis = 1e-9,
      fit_control = list(cv_select = FALSE, prediction_bounds = c(-Inf, Inf),
                         control = list(thresh = 1e-14))
    )
    expect_equal(as.vector(predict(ours, design_matrix(newx, f))),
                 as.vector(predict(theirs, new_data = cbind(newx, -newx))),
                 tolerance = 1e-6, label = paste("symmetric =", sym))
  }
})

# Reflecting one coordinate maps the symmetric dictionary onto itself (up to
# complements, which the free intercept absorbs), so a fit at fixed lambda
# must be the reflection of the original fit. The anchored fit is not.
test_that("symmetric HAL is reflection-equivariant; anchored HAL is not", {
  dat <- and_upper_data()
  xr <- dat$x
  xr[, 2] <- -xr[, 2]
  newx <- toy_x(n = 200, d = ncol(dat$x), seed = 11)
  newxr <- newx
  newxr[, 2] <- -newxr[, 2]
  fit_pred <- function(x, nx, symmetric) {
    f <- make_dictionary(x, 5, symmetric)
    g <- glmnet::glmnet(design_matrix(x, f), dat$y, lambda = 0.01,
                        standardize = FALSE, control = list(thresh = 1e-14))
    as.vector(predict(g, design_matrix(nx, f)))
  }
  expect_equal(fit_pred(dat$x, newx, TRUE), fit_pred(xr, newxr, TRUE),
               tolerance = 1e-5)
  expect_gt(max(abs(fit_pred(dat$x, newx, FALSE) -
                      fit_pred(xr, newxr, FALSE))), 0.05)
})

# With discrete covariates knots fall on data values, where 1{x <= u} is not
# a complement of any "+" main; equivariance needs the "-" mains. Knots are
# spaced wider than the support so the complement's knot is missing.
test_that("symmetric HAL stays reflection-equivariant with tied covariates", {
  set.seed(21)
  x <- matrix(sample(0:20, 600, TRUE), 200, 3)
  y <- 1.5 * (x[, 1] <= 10) + (x[, 2] >= 15) + rnorm(200, sd = 0.3)
  expect_true(any(make_knots(x, 3)[[1]] %in% x[, 1]))
  xr <- x
  xr[, 1] <- -xr[, 1]
  newx <- matrix(sample(0:20, 300, TRUE), 100, 3)
  newxr <- newx
  newxr[, 1] <- -newxr[, 1]
  fit_pred <- function(xx, nx) {
    f <- make_dictionary(xx, 3, TRUE)
    g <- glmnet::glmnet(design_matrix(xx, f), y, lambda = 0.01,
                        standardize = FALSE, control = list(thresh = 1e-14))
    as.vector(predict(g, design_matrix(nx, f)))
  }
  expect_equal(fit_pred(x, newx), fit_pred(xr, newxr), tolerance = 1e-5)
})

# A split sends x >= c one way and x < c the other, so negating a covariate
# swaps the children and leaves the partition, and the fit, unchanged.
test_that("CART is reflection-invariant", {
  dat <- and_upper_data()
  xr <- dat$x
  xr[, 1] <- -xr[, 1]
  grow <- function(x) {
    rpart::rpart(y ~ ., data = data.frame(y = dat$y, x),
                 control = rpart::rpart.control(cp = 0, xval = 0))
  }
  tree <- grow(dat$x)
  expect_gt(sum(tree$frame$var == "x1"), 0)
  expect_equal(unname(predict(tree)), unname(predict(grow(xr))))
})

# Integer covariates put many values exactly halfway between two knots; the
# tie rule has to commute with negation or the first-order symmetric basis
# changes with the coding.
test_that("symmetric knots commute with negation on tied data", {
  source(file.path(root, "R", "nuisance.R"), local = TRUE)
  set.seed(3)
  x <- cbind(age = sample(16:55, 400, TRUE), educ = sample(0:18, 400, TRUE),
             earn = round(rexp(400, 1 / 5000)))
  for (bins in c(50, 25, 12)) {
    xq <- symmetric_knots(x, bins)
    expect_equal(symmetric_knots(-x, bins), -xq)
    expect_true(all(apply(xq, 2, \(v) length(unique(v))) <= bins + 1))
  }
})
