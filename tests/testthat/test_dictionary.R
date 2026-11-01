test_that("design columns equal hand-computed oriented indicator products", {
  x <- toy_x()
  f <- make_dictionary(x, num_knots = 4, symmetric = TRUE)
  design <- as.matrix(design_matrix(x, f))
  by_hand <- vapply(seq_len(nrow(f)), function(i) {
    side <- function(col, u, o) if (o > 0) x[, col] >= u else x[, col] <= u
    v <- side(f$j[i], f$u[i], f$oj[i])
    if (!is.na(f$l[i])) v <- v & side(f$l[i], f$v[i], f$ol[i])
    as.numeric(v)
  }, numeric(nrow(x)))
  expect_equal(unname(design), unname(by_hand))
})

test_that("functions with identical training columns remain available", {
  x <- cbind(1:21, 1:21)
  f <- make_dictionary(x, 3, TRUE)
  expect_equal(nrow(f), 2 * 2 * 3 + 4 * 3^2)
  training <- as.matrix(design_matrix(x, f))
  heldout <- as.matrix(design_matrix(rbind(c(1, 21), c(21, 1)), f))
  first <- which(is.na(f$l) & f$j == 1 & f$u == 11 & f$oj == 1)
  second <- which(is.na(f$l) & f$j == 2 & f$u == 11 & f$oj == 1)
  expect_equal(training[, first], training[, second])
  expect_false(identical(heldout[, first], heldout[, second]))
})

# Columns are compared up to complement: with a free intercept 1 - phi and
# phi span the same model, and a reflected "+" main effect is a "-" one.
test_that("symmetric dictionary is reflection-closed; anchored is not", {
  x <- toy_x(n = 200, d = 3, seed = 7)
  xr <- x
  xr[, 2] <- -xr[, 2]
  canonical <- function(design) {
    m <- as.matrix(design)
    keys <- apply(m, 2, function(col) {
      if (col[1] == 1) col <- 1 - col
      paste(which(col == 1), collapse = ",")
    })
    sort(unique(keys))
  }
  same_span_keys <- function(sym) {
    identical(canonical(design_matrix(x, make_dictionary(x, 5, sym))),
              canonical(design_matrix(xr, make_dictionary(xr, 5, sym))))
  }
  expect_true(same_span_keys(TRUE))
  expect_false(same_span_keys(FALSE))
})

test_that("the anchored dictionary is a subset of the symmetric one", {
  x <- toy_x()
  a <- make_dictionary(x, 4, FALSE)
  s <- make_dictionary(x, 4, TRUE)
  key <- function(f) paste(f$j, f$l, f$u, f$v, f$oj, f$ol)
  expect_true(all(key(a) %in% key(s)))
})

test_that("closed discrete indicators can enlarge the prediction space", {
  x <- as.matrix(expand.grid(x1 = (0:4) / 4, x2 = (0:4) / 4))
  knots <- rep(list(c(0.25, 0.5, 0.75)), 2)
  design <- function(x, knots, symmetric) {
    cbind(1, as.matrix(design_matrix(
      x, dictionary_features(knots, symmetric)
    )))
  }
  anchored <- design(x, knots, FALSE)
  symmetric <- design(x, knots, TRUE)
  expect_equal(qr(anchored)$rank, 16)
  expect_equal(qr(symmetric)$rank, 25)
  expect_equal(anchored[x[, 1] == 0.75, ], anchored[x[, 1] == 1, ])
  expect_false(identical(symmetric[x[, 1] == 0.75, ],
                         symmetric[x[, 1] == 1, ]))
  duplicate <- cbind(x, x[, 2])
  knots <- rep(knots[1], 3)
  expect_equal(qr(design(duplicate, knots, FALSE))$rank, 16)
  duplicate[, 3] <- -duplicate[, 3]
  knots[[3]] <- -knots[[3]]
  expect_equal(qr(design(duplicate, knots, FALSE))$rank, 20)
  expect_equal(qr(design(duplicate, knots, TRUE))$rank, 25)
})
