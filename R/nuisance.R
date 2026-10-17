# HAL fitters shared by the real-data experiments. Each fits on (x, y) with
# the given CV folds and returns predictions (on the response scale) at each
# matrix in the list `newx`.

hal_cv <- function(x, y, symmetric, family, foldid, newx) {
  features <- make_dictionary(x, 9, symmetric)
  # The default convergence threshold leaves CV curves loose enough that
  # lambda.min can differ between equivalent codings; tighten it so any
  # remaining difference is the estimator's, not the solver's.
  cv <- glmnet::cv.glmnet(design_matrix(x, features), y, family = family,
                          standardize = FALSE, foldid = foldid,
                          control = list(thresh = 1e-12))
  map(newx, \(nx) {
    as.vector(predict(cv, design_matrix(nx, features), s = "lambda.min",
                      type = "response"))
  })
}

# hal9001 fits keep glmnet's default tolerance: at 1e-12 glmnet stops before
# the smallest lambdas on these designs, which truncates the CV path.
hal9001_cv <- function(x, y, family, foldid, newx) {
  fit <- hal9001::fit_hal(x, y, family = family,
                          fit_control = list(foldid = foldid))
  map(newx, \(nx) as.vector(predict(fit, new_data = nx)))
}

# Symmetric first-order HAL through hal9001: the design is cbind(x, -x), so
# hinges on -x_j are the reflected hinges (u - x_j)_+, and terms that use both
# copies of a covariate are dropped. hal9001's own knot rule keeps the minimum
# and drops the top 2% of each covariate, which is not reflection symmetric, so
# knots here are type-7 quantiles at levels evenly spaced on [0, 1] with each
# value snapped to its nearest knot. Knot counts per degree are hal9001's.
symmetric_knots <- function(x, bins) {
  apply(x, 2, \(v) {
    if (length(unique(v)) <= bins) return(v)
    q <- unique(stats::quantile(v, seq(0, 1, length.out = bins), type = 7))
    dist <- abs(outer(v, q, `-`))
    # A value halfway between two knots goes to the one nearer the median:
    # negation swaps the two knots and negates the median, so the rule
    # commutes with reflection, where "the lower knot" would not. A value
    # that is itself the midpoint of two knots around the median stays put.
    tol <- 1e-9 * diff(range(v))
    near <- dist <= apply(dist, 1, min) + tol
    score <- ifelse(near, abs(rep(q, each = length(v)) - stats::median(v)),
                    Inf)
    best <- score <= apply(score, 1, min) + tol
    ifelse(rowSums(best) > 1, v, q[max.col(-score, "first")])
  })
}

hal9001_sym_cv <- function(x, y, family, foldid, newx) {
  k <- ncol(x)
  num_knots <- c(50, 25, 12)
  basis <- unlist(map(1:3, \(degree) {
    xq <- symmetric_knots(x, num_knots[degree])
    b <- hal9001::enumerate_basis(cbind(xq, -xq), max_degree = degree,
                                  smoothness_orders = rep(1, 2 * k))
    keep(b, \(e) length(e$cols) == degree && !anyDuplicated((e$cols - 1) %% k))
  }), recursive = FALSE)
  fit <- hal9001::fit_hal(cbind(x, -x), y, family = family,
                          basis_list = unique(basis), smoothness_orders = 1,
                          max_degree = 3, fit_control = list(foldid = foldid))
  map(newx, \(nx) as.vector(predict(fit, new_data = cbind(nx, -nx))))
}

# First-order hinges carry the units of their covariate and hal9001 turns
# off glmnet's standardization, so the penalty depends on units. Dividing
# each covariate by its training standard deviation removes that; shifts
# already leave the fit unchanged.
hal9001_sym_scaled_cv <- function(x, y, family, foldid, newx) {
  spread <- apply(x, 2, stats::sd)
  spread[spread == 0] <- 1
  hal9001_sym_cv(sweep(x, 2, spread, `/`), y, family, foldid,
                 map(newx, \(nx) sweep(nx, 2, spread, `/`)))
}

fit_nuisance <- function(method, x, y, family, foldid, newx) {
  switch(method,
    hal = hal_cv(x, y, FALSE, family, foldid, newx),
    shal = hal_cv(x, y, TRUE, family, foldid, newx),
    hal9001 = hal9001_cv(x, y, family, foldid, newx),
    shal9001 = hal9001_sym_cv(x, y, family, foldid, newx),
    shal9001_scaled = hal9001_sym_scaled_cv(x, y, family, foldid, newx)
  )
}
