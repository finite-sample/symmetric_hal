test_that("nominal label permutations change a saturated threshold penalty", {
  x <- matrix(0:2, ncol = 1)
  recoded <- matrix(c(1, 0, 2), ncol = 1)
  y <- c(0, 1, 0)
  first <- fit_hal(x, y, 1 / 18, knots = list(0:2))
  second <- fit_hal(recoded, y, 1 / 18, knots = list(0:2))
  expect_equal(qr(cbind(1, design_matrix(x, first$features)))$rank, 3)
  expect_equal(qr(cbind(1, design_matrix(recoded, second$features)))$rank, 3)
  expect_equal(first$fitted, c(1 / 6, 2 / 3, 1 / 6))
  expect_equal(second$fitted, c(1 / 12, 5 / 6, 1 / 12))
  expect_equal(sum(abs(first$beta)), 1)
  expect_equal(sum(abs(second$beta)), 3 / 4)
  difference <- rbind(c(-1, 1, 0), c(0, -1, 1))
  first_score <- (first$fitted - y) / 3 +
    crossprod(difference, c(1, -1)) / 18
  expect_equal(as.vector(first_score), rep(0, 3))
  ordering <- order(recoded[, 1])
  second_score <- (second$fitted[ordering] - y[ordering]) / 3 +
    crossprod(difference, c(-1, -1 / 2)) / 18
  expect_equal(as.vector(second_score), rep(0, 3))
})

test_that("full category indicators preserve tuned fits under relabeling", {
  x <- diag(4)[rep_len(1:3, 60), ]
  newx <- diag(4)
  permutation <- c(3, 1, 4, 2)
  folds <- rep_len(1:5, nrow(x))
  responses <- list(
    gaussian = x[, 2] + (seq_len(nrow(x)) %% 7) / 10,
    binomial = as.numeric(seq_len(nrow(x)) %% 4 < 2)
  )
  for (symmetric in c(FALSE, TRUE)) {
    for (family in names(responses)) {
      first <- cv_hal(x, responses[[family]], c(0.05, 0.2), folds,
                      family = family, symmetric = symmetric, num_knots = 1)
      second <- cv_hal(x[, permutation], responses[[family]],
                       c(0.05, 0.2), folds, family = family,
                       symmetric = symmetric, num_knots = 1)
      expect_equal(first$lambda, second$lambda)
      expect_equal(first$cv$loss, second$cv$loss, tolerance = 1e-7)
      expect_equal(first$predict(newx), second$predict(newx[, permutation]),
                   tolerance = 1e-6)
    }
  }
})

test_that("zero-order tuned fits respect affine covariate recodings", {
  x <- as.matrix(expand.grid(x1 = (0:6) / 4, x2 = -3:3))
  newx <- rbind(c(0.1, -1.3), c(0.8, 0.2), c(1.2, 2.5))
  transform <- function(z) {
    sweep(sweep(z, 2, c(8, -0.25), `*`), 2, c(16, 2), `+`)
  }
  folds <- rep_len(1:5, nrow(x))
  responses <- list(
    gaussian = x[, 1] + (x[, 2] >= 0) + (seq_len(nrow(x)) %% 3) / 10,
    binomial = as.numeric(seq_len(nrow(x)) %% 3 == 0)
  )
  for (family in names(responses)) {
    fit <- cv_hal(x, responses[[family]], c(0.05, 0.2), folds,
                  family = family)
    other <- cv_hal(transform(x), responses[[family]], c(0.05, 0.2), folds,
                    family = family)
    expect_equal(other$lambda, fit$lambda)
    expect_equal(other$cv$loss, fit$cv$loss, tolerance = 1e-7)
    expect_equal(other$predict(transform(newx)), fit$predict(newx),
                 tolerance = 1e-6)
  }
})

test_that("interpolated knots have the nonlinear-recoding counterexample", {
  x <- matrix(c(0, 2), ncol = 1)
  y <- c(0, 1)
  features <- make_dictionary(x, 3, TRUE)
  recoded <- make_dictionary(x^3, 3, TRUE)
  design <- as.matrix(design_matrix(x, features))
  other_design <- as.matrix(design_matrix(x^3, recoded))
  expect_equal(design, other_design)
  beta <- ifelse(features$oj > 0, 2 / 15, -2 / 15)
  fitted <- 0.5 + as.vector(design %*% beta)
  gradient <- as.vector(crossprod(design, fitted - y)) / length(y)
  expect_equal(fitted, c(0.1, 0.9))
  expect_equal(gradient + 0.05 * sign(beta), rep(0, length(beta)),
               tolerance = 1e-14)
  expect_equal(mean(fitted - y), 0, tolerance = 1e-14)
  expect_equal(sum(beta^2), (4 / 5)^2 / 6)
  original <- 0.5 + as.vector(design_matrix(matrix(3 / 4), features) %*% beta)
  transformed <- 0.5 + as.vector(
    design_matrix(matrix((3 / 4)^3), recoded) %*% beta
  )
  expect_equal(original, 11 / 30)
  expect_equal(transformed, 1 / 10)
})

test_that("rank knots preserve functions under monotone recodings", {
  samples <- list(
    cbind(c(-2, -1, 0, 1, 2, 3), c(1, 3, 2, -1, 0, -2)),
    cbind(c(-2, -2, 0, 0, 2, 2), c(1, 1, 1, 1, 1, 1))
  )
  recode <- function(z) cbind(exp(z[, 1]), -z[, 2]^3)
  rank_knots <- function(z) {
    lapply(seq_len(ncol(z)), function(j) {
      unique(sort(z[, j])[c(2, 3, 4, 5)])
    })
  }
  feature_key <- function(features) do.call(paste, c(features, sep = "|"))
  queries <- as.matrix(expand.grid(
    c(-3, -2, -1.5, 0, 0.5, 2, 3, 4),
    c(-3, -2, -0.5, 0, 1, 1.5, 3, 4)
  ))
  for (x in samples) {
    knots <- rank_knots(x)
    other_knots <- rank_knots(recode(x))
    features <- dictionary_features(knots, TRUE)
    other_features <- dictionary_features(other_knots, TRUE)
    mapped <- features
    mapped$u <- ifelse(mapped$j == 1, exp(mapped$u), -mapped$u^3)
    mapped$oj <- mapped$oj * ifelse(mapped$j == 1, 1, -1)
    interaction <- !is.na(mapped$l)
    mapped$v[interaction] <- -mapped$v[interaction]^3
    mapped$ol[interaction] <- -mapped$ol[interaction]
    permutation <- match(feature_key(mapped), feature_key(other_features))
    expect_equal(sort(permutation), seq_len(nrow(features)))
    expect_equal(
      as.matrix(design_matrix(queries, features)),
      as.matrix(design_matrix(recode(queries), other_features))[, permutation]
    )
    fit_rows <- rep(seq_len(nrow(x)), each = 3)
    fit_x <- x[fit_rows, , drop = FALSE]
    for (family in c("gaussian", "binomial")) {
      y <- if (family == "gaussian") {
        c(0, 1, 2, 1, -1, 0)
      } else {
        c(0, 1, 0, 1, 1, 0)
      }
      fit <- fit_hal(fit_x, y[fit_rows], 0.05, family = family, knots = knots)
      other <- fit_hal(recode(fit_x), y[fit_rows], 0.05, family = family,
                       knots = other_knots)
      expect_equal(other$predict(recode(queries)), fit$predict(queries),
                   tolerance = 1e-6)
    }
  }
})

test_that("observed knots keep both central values and commute with reversal", {
  expect_equal(make_knots(matrix(1:6), 1, "observed"), list(c(3L, 4L)))
  expect_equal(make_knots(matrix(1:5), 1, "observed"), list(3L))
  expect_equal(make_knots(matrix(c(0, 2)), 3, "observed"), list(c(0, 2)))
  expect_equal(make_knots(matrix(rep(7, 6)), 3, "observed"), list(7))
  for (n in 2:12) {
    for (q in 1:7) {
      x <- matrix(seq_len(n))
      knots <- make_knots(x, q, "observed")[[1]]
      expect_equal(sort(n + 1 - knots), knots)
      expect_true(all(knots %in% x))
      expect_lte(length(knots), 2 * q)
    }
  }
})

test_that("observed knots give fold-local nonlinear coding equivariance", {
  x <- as.matrix(expand.grid(x1 = (0:6) / 4, x2 = -3:3))
  newx <- rbind(c(0.1, -1.3), c(0.8, 0.2), c(1.2, 2.5), c(-1, 4))
  recode <- function(z) cbind(exp(z[, 1]), -z[, 2]^3)
  folds <- rep_len(1:5, nrow(x))
  responses <- list(
    gaussian = x[, 1] + (x[, 2] >= 0) + (seq_len(nrow(x)) %% 3) / 10,
    binomial = as.numeric(seq_len(nrow(x)) %% 3 == 0)
  )
  for (family in names(responses)) {
    y <- responses[[family]]
    fit <- cv_hal(x, y, c(0.05, 0.2), folds,
                  family = family, knot_rule = "observed")
    other <- cv_hal(recode(x), y, c(0.05, 0.2), folds,
                    family = family, knot_rule = "observed")
    expect_equal(fit$knot_rule, "observed")
    expect_equal(other$lambda, fit$lambda)
    expect_equal(other$cv$loss, fit$cv$loss, tolerance = 1e-7)
    expect_equal(other$validation_predictions, fit$validation_predictions,
                 tolerance = 1e-6)
    expect_equal(other$predict(recode(newx)), fit$predict(newx),
                 tolerance = 1e-6)
    for (fold in unique(folds)) {
      validation <- folds == fold
      manual <- fit_hal(x[!validation, ], y[!validation], fit$cv$lambda[1],
                        family = family, knot_rule = "observed")
      expect_equal(fit$validation_predictions[validation, 1],
                   manual$predict(x[validation, ], type = "link"))
    }
  }
  expect_error(fit_hal(x, responses$gaussian, 0.1, knot_rule = "unknown"),
               "arg")
  expect_error(fit_hal(x, responses$gaussian, 0.1,
                       knots = list(0.5, 0), knot_rule = "observed"),
               "supplied knots")
})
