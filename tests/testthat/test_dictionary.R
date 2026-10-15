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

test_that("no constant or duplicate columns", {
  x <- toy_x()
  for (sym in c(FALSE, TRUE)) {
    f <- make_dictionary(x, 4, sym)
    design <- design_matrix(x, f)
    share <- Matrix::colMeans(design)
    expect_true(all(share > 0 & share < 1))
    expect_false(anyDuplicated(column_keys(design)) > 0)
  }
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
