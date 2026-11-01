test_that("anchor selection minimizes joint validation loss", {
  source(file.path(root, "R", "anchors.R"), local = TRUE)
  x <- toy_x(32, 2, 241)
  y <- as.numeric(x[, 1] < 0 & x[, 2] > 0) + sin(1:32) / 5
  folds <- rep(1:4, 8)
  fit <- cv_anchors(x, y, c(0.02, 0.1), folds, num_knots = 2)
  expect_equal(nrow(fit$signs), 4)
  selected_loss <- with(fit$cv, loss[anchor == fit$selected_anchor &
                                       lambda == fit$selected_lambda])
  expect_lte(selected_loss, min(fit$cv$loss) + 1e-8)
  newx <- toy_x(12, 2, 242)
  separate <- cv_hal(x, y, c(0.02, 0.1), folds,
                     num_knots = 2, symmetric = FALSE)
  expect_equal(fit$predict_each(newx)[, 1], separate$predict(newx))
  for (j in seq_along(fit$average_cv$lambda)) {
    prediction <- Reduce(`+`, lapply(fit$validation_predictions, function(p) {
      p[, j]
    })) / 4
    expect_equal(fit$average_cv$loss[j], mean((y - prediction)^2) / 2)
  }
  expect_true(all(fit$seconds >= 0))
})

test_that("anchor averaging is invariant under reflection with shared folds", {
  source(file.path(root, "R", "anchors.R"), local = TRUE)
  x <- cbind(1:24, 1:24)
  folds <- rep(1:3, 8)
  newx <- rbind(c(1, 24), c(24, 1), c(10, 20))
  for (family in c("gaussian", "binomial")) {
    y <- as.numeric(x[, 1] > 12)
    original <- cv_anchors(x, y, c(0.03, 0.1), folds,
                           num_knots = 2, family = family)
    reflected <- cv_anchors(sweep(x, 2, c(-1, 1), `*`), y,
                            c(0.03, 0.1), folds,
                            num_knots = 2, family = family)
    expect_equal(original$average_cv$loss, reflected$average_cv$loss,
                 tolerance = 1e-7)
    expect_equal(original$predict_average(newx),
                 reflected$predict_average(sweep(newx, 2, c(-1, 1), `*`)),
                 tolerance = 1e-7)
    average <- rowMeans(vapply(seq_len(4), function(i) {
      original$average_fits[[i]]$predict(sweep(newx, 2,
                                               original$signs[i, ], `*`))
    }, numeric(nrow(newx))))
    expect_equal(original$predict_average(newx), average)
    expect_equal(original$predict_average(newx[1, , drop = FALSE]), average[1])
  }
})

test_that("flat validation losses select largest penalty then first anchor", {
  source(file.path(root, "R", "anchors.R"), local = TRUE)
  fit <- cv_anchors(toy_x(24, 2, 243), rep(3, 24), c(0.1, 1), rep(1:3, 8))
  expect_equal(fit$selected_lambda, 1)
  expect_equal(fit$average_lambda, 1)
  expect_equal(fit$selected_anchor, 1)
})
