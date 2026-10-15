# Independent evaluator for features of any degree, used to check the
# propositions without going through R/basis.R.
eval_feature <- function(x, cols, knots, orient) {
  v <- rep(TRUE, nrow(x))
  for (m in seq_along(cols)) {
    v <- v & if (orient[m] > 0) x[, cols[m]] >= knots[m] else
      x[, cols[m]] <= knots[m]
  }
  as.numeric(v)
}

full_dictionary <- function(x, grid, symmetric, max_degree) {
  d <- ncol(x)
  out <- list()
  for (r in seq_len(max_degree)) {
    for (s in utils::combn(d, r, simplify = FALSE)) {
      knot_sets <- as.matrix(expand.grid(rep(list(grid), r)))
      orients <- matrix(1, 1, r)
      if (symmetric) orients <- as.matrix(expand.grid(rep(list(c(1, -1)), r)))
      for (a in seq_len(nrow(knot_sets))) for (o in seq_len(nrow(orients))) {
        out[[length(out) + 1]] <- eval_feature(x, s, knot_sets[a, ],
                                               orients[o, ])
      }
    }
  }
  do.call(cbind, out)
}

# min ||beta||_1 subject to b0 + F beta = target, as an LP in
# (beta+, beta-, b0+, b0-) >= 0, solved exactly by the simplex method.
min_l1 <- function(features, target) {
  p <- ncol(features)
  a <- cbind(features, -features, 1, -1)
  sol <- lpSolve::lp("min", c(rep(1, 2 * p), 0, 0), a,
                     rep("=", nrow(a)), target)
  stopifnot(sol$status == 0)
  sol$objval
}

# Symmetric: 2^t. Anchored: 2^(t+k) - 1{m = t = 0}.
test_that("Prop 1: box penalties match the closed forms", {
  grid <- c(0.2, 0.4, 0.6, 0.8)
  mids <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- as.matrix(expand.grid(mids, mids, mids))
  a_dict <- full_dictionary(x, grid, FALSE, 3)
  s_dict <- full_dictionary(x, grid, TRUE, 3)
  set.seed(11)
  for (trial in 1:12) {
    s <- sort(sample(3, sample(1:3, 1)))
    type <- sample(c("upper", "lower", "two"), length(s), replace = TRUE)
    inside <- rep(TRUE, nrow(x))
    for (m in seq_along(s)) {
      ab <- sort(sample(grid, 2))
      v <- x[, s[m]]
      inside <- inside & switch(type[m], upper = v >= ab[1],
                                lower = v <= ab[2],
                                two = v >= ab[1] & v <= ab[2])
    }
    cval <- 1.7
    target <- cval * inside
    t_two <- sum(type == "two")
    k_low <- sum(type == "lower")
    m_up <- sum(type == "upper")
    expected_a <- cval * (2^(t_two + k_low) - (m_up == 0 && t_two == 0))
    expected_s <- cval * 2^t_two
    # The anchored representation is unique, so least squares recovers it.
    coefs <- qr.coef(qr(cbind(1, a_dict)), target)
    expect_equal(sum(abs(coefs[-1]), na.rm = TRUE), expected_a,
                 tolerance = 1e-8, label = paste(type, collapse = "/"))
    expect_equal(min_l1(s_dict, target), expected_s, tolerance = 1e-9,
                 label = paste(type, collapse = "/"))
  }
})

test_that("Thm 1: oracle inequality holds under its lambda condition", {
  set.seed(5)
  x <- toy_x(n = 300, d = 3, seed = 5)
  f <- make_dictionary(x, 5, FALSE)
  design <- as.matrix(design_matrix(x, f))
  fstar <- 1.5 * (x[, 1] >= 0) * (x[, 2] >= 0.2)
  eps <- rnorm(300, sd = 0.7)
  y <- fstar + eps
  pc <- function(v) v - mean(v)
  px <- apply(design, 2, pc)
  lambda <- 2 * max(abs(crossprod(px, eps))) / 300
  fit <- glmnet::glmnet(design, y, lambda = lambda, standardize = FALSE,
                        control = list(thresh = 1e-14))
  fhat <- as.vector(predict(fit, design))
  lhs <- mean(pc(fhat - fstar)^2)
  # Oracle beta: least-squares fit of fstar on the dictionary (exact if
  # representable, best approximation otherwise).
  ls <- lm.fit(cbind(1, design), fstar)
  beta <- ls$coefficients[-1]
  beta[is.na(beta)] <- 0
  rhs <- mean(pc(cbind(1, design) %*% c(ls$coefficients[1], beta) - fstar)^2) +
    3 * lambda * sum(abs(beta))
  expect_lte(lhs, rhs)
})
