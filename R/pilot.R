pilot_truths <- list(
  upper = function(x) as.numeric(x[, 1] >= 0.75 & x[, 2] >= 0.75),
  lower = function(x) as.numeric(x[, 1] <= 0.25 & x[, 2] <= 0.25),
  mixed = function(x) as.numeric(x[, 1] >= 0.75 & x[, 2] <= 0.25),
  interior = function(x) {
    as.numeric(abs(x[, 1] - 0.5) <= 0.25 & abs(x[, 2] - 0.5) <= 0.25)
  },
  conflicting = function(x) {
    as.numeric(x[, 1] <= 0.25 & x[, 2] <= 0.25) +
      as.numeric(x[, 1] >= 0.75 & x[, 2] >= 0.75)
  },
  additive = function(x) {
    as.numeric(x[, 1] >= 0.25) + 2 * as.numeric(x[, 2] <= 0.75)
  },
  smooth = function(x) sin(2 * pi * x[, 1]) + 4 * (x[, 2] - 0.5)^2
)

pilot_x <- function(n, d, distribution) {
  x <- switch(distribution,
    independent = matrix(runif(n * d), n, d),
    correlated = {
      z <- matrix(rnorm(n * d), n, d)
      pnorm(z %*% chol(0.5 + 0.5 * diag(d)))
    },
    discrete = {
      z <- matrix(sample(0:4, n * d, replace = TRUE) / 4, n, d)
      if (d == 3) z[, 3] <- z[, 2]
      z
    },
    stop("unknown covariate distribution")
  )
  colnames(x) <- paste0("x", seq_len(d))
  x
}

pilot_grid <- function(reps = 20) {
  grid <- expand.grid(
    truth = names(pilot_truths), d = c(2, 3), n = c(100, 300),
    snr = c(0.5, 2), distribution = c("independent", "correlated", "discrete"),
    replicate = seq_len(reps), stringsAsFactors = FALSE
  )
  grid$knot_rule <- "quantile"
  fixed <- subset(grid, truth == "conflicting" & d == 2 &
                    distribution == "independent" & snr == 2)
  fixed$knot_rule <- "fixed"
  grid <- rbind(grid, fixed)
  grid <- grid[order(grid$replicate, grid$d, grid$n, grid$snr), ]
  grid$job <- seq_len(nrow(grid))
  grid
}

pilot_moments <- function(x, truth, lower = NULL, upper = NULL) {
  f <- pilot_truths[[truth]](x)
  if (truth != "smooth" || is.null(lower)) {
    return(list(mean = f, second = f^2))
  }
  a <- lower[, 1]
  b <- upper[, 1]
  s1 <- (cos(2 * pi * a) - cos(2 * pi * b)) / (2 * pi * (b - a))
  s2 <- 0.5 - (sin(4 * pi * b) - sin(4 * pi * a)) / (8 * pi * (b - a))
  a <- lower[, 2] - 0.5
  b <- upper[, 2] - 0.5
  m2 <- (b^3 - a^3) / (3 * (b - a))
  m4 <- (b^5 - a^5) / (5 * (b - a))
  list(mean = s1 + 4 * m2, second = s2 + 8 * s1 * m2 + 16 * m4)
}

pilot_integration <- function(x, truth, distribution, knots) {
  d <- ncol(x)
  if (distribution == "correlated") {
    f <- pilot_truths[[truth]](x)
    return(list(x = x, weights = rep(1 / nrow(x), nrow(x)),
                mean = f, second = f^2, exact = FALSE))
  }
  if (distribution == "discrete") {
    grid <- as.matrix(expand.grid(rep(list((0:4) / 4), 2)))
    if (d == 3) grid <- cbind(grid, grid[, 2])
    f <- pilot_truths[[truth]](grid)
    return(list(x = grid, weights = rep(1 / nrow(grid), nrow(grid)),
                mean = f, second = f^2, exact = TRUE))
  }
  cuts <- lapply(knots, function(k) sort(unique(c(0, k, 0.25, 0.75, 1))))
  intervals <- lapply(cuts, function(k) seq_len(length(k) - 1))
  indices <- as.matrix(expand.grid(intervals))
  lower <- upper <- indices * 0
  for (j in seq_len(d)) {
    lower[, j] <- cuts[[j]][indices[, j]]
    upper[, j] <- cuts[[j]][indices[, j] + 1]
  }
  grid <- (lower + upper) / 2
  weights <- apply(upper - lower, 1, prod)
  moments <- pilot_moments(grid, truth, lower, upper)
  list(x = grid, weights = weights, mean = moments$mean,
       second = moments$second, exact = TRUE)
}

pilot_draw <- function(cell, signal_variance) {
  truth_id <- match(cell$truth, names(pilot_truths))
  set.seed(50000 + cell$replicate * 1000 + truth_id * 10 +
             cell$d + cell$n + as.integer(cell$snr * 100))
  x <- pilot_x(cell$n, cell$d, cell$distribution)
  y <- pilot_truths[[cell$truth]](x) +
    rnorm(cell$n, sd = sqrt(signal_variance / cell$snr))
  folds <- sample(rep_len(1:5, cell$n))
  list(x = x, y = y, folds = folds)
}

pilot_replicate <- function(cell, test_x, signal_variance,
                            lambda = c(0.003, 0.01, 0.03, 0.1, 0.3)) {
  started <- proc.time()[["elapsed"]]
  draw <- pilot_draw(cell, signal_variance)
  x <- draw$x
  y <- draw$y
  folds <- draw$folds
  centre <- mean(y)
  spread <- sd(y)
  z <- (y - centre) / spread
  fixed <- if (cell$knot_rule == "fixed") rep(list(c(0.25, 0.5, 0.75)), cell$d)
  suite <- cv_anchors(x, z, lambda, folds, knots = fixed)
  symmetric <- cv_hal(x, z, lambda, folds, knots = fixed)
  grid_knots <- if (is.null(fixed)) make_knots(x, 3) else fixed
  integration <- pilot_integration(test_x, cell$truth, cell$distribution,
                                   grid_knots)
  predict_all <- function(newx) {
    each <- suite$predict_each(newx)
    cbind(anchored = each[, 1], selected = suite$predict_selected(newx),
          averaged = suite$predict_average(newx),
          symmetric = symmetric$predict(newx), each)
  }
  prediction <- centre + spread * predict_all(integration$x)
  risk <- colSums(integration$weights *
                    (prediction^2 - 2 * prediction * integration$mean +
                       integration$second))
  oracle <- which.min(risk[-(1:4)])
  methods <- c("anchored", "selected", "averaged", "symmetric", "oracle_anchor")
  test_prediction <- centre + spread * predict_all(test_x)
  test_truth <- pilot_truths[[cell$truth]](test_x)
  test_risk <- colMeans((test_prediction - test_truth)^2)
  diagnostics <- do.call(rbind, lapply(seq_along(suite$fits), function(i) {
    fit <- suite$fits[[i]]
    rbind(data.frame(anchor = i, fit$cv_diagnostics),
          data.frame(anchor = i, fold = 0, lambda = fit$lambda,
                     t(fit$diagnostics), check.names = FALSE))
  }))
  sym_diagnostics <- rbind(symmetric$cv_diagnostics,
                           data.frame(fold = 0, lambda = symmetric$lambda,
                                      t(symmetric$diagnostics)))
  signs <- anchor_signs(cell$d)
  movement <- numeric(nrow(signs))
  reflection_fits <- list()
  original_test <- symmetric$predict(test_x)
  for (i in seq_len(nrow(signs))[-1]) {
    reflected_knots <- if (!is.null(fixed)) {
      lapply(seq_len(cell$d), function(j) signs[i, j] * fixed[[j]])
    }
    fit <- cv_hal(sweep(x, 2, signs[i, ], `*`), z, lambda, folds,
                  knots = reflected_knots)
    movement[i] <- max(abs(fit$predict(sweep(test_x, 2, signs[i, ], `*`)) -
                             original_test)) * spread
    reflection_fits[[i]] <- list(lambda = fit$lambda, cv = fit$cv,
                                 cv_diagnostics = fit$cv_diagnostics,
                                 diagnostics = fit$diagnostics)
  }
  if (max(movement) > 1e-6 * max(1, spread)) {
    stop("held-out reflection check failed: ", max(movement))
  }
  rows <- data.frame(
    method = methods, risk = c(risk[1:4], risk[4 + oracle]),
    test_risk = c(test_risk[1:4], test_risk[4 + oracle]),
    lambda = c(suite$fits[[1]]$lambda, suite$selected_lambda,
               suite$average_lambda, symmetric$lambda,
               suite$fits[[oracle]]$lambda),
    seconds = c(suite$seconds, sum(symmetric$seconds), NA_real_),
    signal_variance = signal_variance, exact_integration = integration$exact,
    selected_anchor = suite$selected_anchor,
    max_symmetric_movement = max(movement), outcome_sd = spread
  )
  for (method in c("cart", "rf")) {
    start <- proc.time()[["elapsed"]]
    set.seed(90000 + cell$replicate)
    fit <- if (method == "cart") fit_cart(x, y) else fit_rf(x, y)
    prediction <- fit$predict(test_x)
    elapsed <- proc.time()[["elapsed"]] - start
    row <- rows[1, ]
    row$method <- method
    row$risk <- row$test_risk <- mean((prediction - test_truth)^2)
    row$lambda <- NA_real_
    row$seconds <- elapsed
    row$exact_integration <- FALSE
    rows <- rbind(rows, row)
  }
  stopifnot(all(is.finite(rows$risk)), all(rows$risk >= -1e-10))
  list(cell = cell, results = rows,
       elapsed_seconds = proc.time()[["elapsed"]] - started,
       anchor_cv = suite$cv,
       average_cv = suite$average_cv, symmetric_cv = symmetric$cv,
       anchor_diagnostics = diagnostics,
       symmetric_diagnostics = sym_diagnostics,
       reflection_movement = movement, reflection_fits = reflection_fits,
       average_refit_diagnostics = lapply(suite$average_fits, function(fit) {
         fit$diagnostics
       }),
       selected_refit_diagnostics = suite$selected_fit$diagnostics,
       foldid = folds)
}
