and_upper_data <- function(n = 200, seed = 4) {
  x <- toy_x(n = n, d = 3, seed = seed)
  set.seed(seed + 1)
  list(x = x, y = 1.5 * (x[, 1] >= 0.2 & x[, 2] >= 0.2) + rnorm(n, sd = 0.3))
}

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
