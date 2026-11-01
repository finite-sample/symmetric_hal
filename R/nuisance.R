# Shared nuisance fits use the same explicit folds and candidate penalties.
hal_cv <- function(x, y, symmetric, family, foldid, newx, canonical = FALSE) {
  centre <- if (family == "gaussian") mean(y) else 0
  spread <- if (family == "gaussian") max(stats::sd(y), 1e-8) else 1
  response <- (y - centre) / spread
  lambda <- if (family == "gaussian") {
    c(0.003, 0.01, 0.03, 0.1)
  } else {
    c(0.002, 0.01, 0.05, 0.2)
  }
  fit <- if (canonical) {
    cv_canonical_hal(x, response, lambda, foldid, family = family)
  } else {
    cv_hal(x, response, lambda, foldid, symmetric = symmetric, family = family)
  }
  result <- lapply(newx, function(z) centre + spread * fit$predict(z))
  attr(result, "diagnostics") <- list(
    lambda = fit$lambda, cv = fit$cv, certificate = fit$diagnostics,
    cv_certificates = fit$cv_diagnostics, outcome_centre = centre,
    outcome_scale = spread,
    basis_size = if (canonical) {
      sum(vapply(fit$fits, function(part) nrow(part$features), integer(1)))
    } else {
      nrow(fit$features)
    },
    orientations = if (canonical) fit$signs else NULL
  )
  result
}

hal9001_cv <- function(x, y, family, foldid, newx, scaled = FALSE,
                       tight = FALSE) {
  if (scaled) {
    centre <- colMeans(x)
    spread <- apply(x, 2, stats::sd)
    spread[spread == 0] <- 1
    transform <- function(z) sweep(sweep(z, 2, centre, `-`), 2, spread, `/`)
    x <- transform(x)
    newx <- lapply(newx, transform)
  }
  control <- list(foldid = foldid)
  if (tight) control$control <- list(maxit = 1000000, thresh = 1e-12)
  fit <- hal9001::fit_hal(x, y, family = family,
                          fit_control = control)
  result <- lapply(newx, \(nx) as.vector(predict(fit, new_data = nx)))
  attr(result, "diagnostics") <- list(lambda = fit$lambda_star,
                                      cv = fit$lasso_fit$cvm,
                                      lambda_grid = fit$lasso_fit$lambda,
                                      basis_size = length(fit$basis_list))
  result
}

fit_nuisance <- function(method, x, y, family, foldid, newx) {
  switch(method,
    hal = hal_cv(x, y, FALSE, family, foldid, newx),
    shal = hal_cv(x, y, TRUE, family, foldid, newx),
    canonical = hal_cv(x, y, FALSE, family, foldid, newx, canonical = TRUE),
    hal9001 = hal9001_cv(x, y, family, foldid, newx),
    hal9001_canonical = hal9001_canonical_cv(x, y, family, foldid, newx),
    hal9001_tight = hal9001_cv(x, y, family, foldid, newx, tight = TRUE),
    hal9001_scaled = hal9001_cv(x, y, family, foldid, newx, scaled = TRUE),
    hal9001_scaled_tight = hal9001_cv(x, y, family, foldid, newx,
                                      scaled = TRUE, tight = TRUE),
    stop("unknown nuisance method: ", method)
  )
}

hal9001_canonical_cv <- function(x, y, family, foldid, newx) {
  pivot <- apply(x, 2, stats::median)
  residual <- sweep(x, 2, pivot, `-`)
  centre <- colMeans(residual)
  spread <- apply(residual, 2, stats::sd)
  spread[spread == 0] <- 1
  transform <- function(z) {
    sweep(sweep(sweep(z, 2, pivot, `-`), 2, centre, `-`), 2, spread, `/`)
  }
  x <- transform(x)
  newx <- lapply(newx, transform)
  signs <- canonical_signs(x)
  fits <- lapply(seq_len(nrow(signs)), function(i) {
    hal9001_cv(sweep(x, 2, signs[i, ], `*`), y, family, foldid,
               lapply(newx, function(z) sweep(z, 2, signs[i, ], `*`)))
  })
  result <- lapply(seq_along(newx), function(j) {
    Reduce(`+`, lapply(fits, `[[`, j)) / length(fits)
  })
  diagnostics <- lapply(fits, attr, which = "diagnostics")
  attr(result, "diagnostics") <- list(
    lambda = vapply(diagnostics, `[[`, numeric(1), "lambda"),
    basis_size = vapply(diagnostics, `[[`, integer(1), "basis_size"),
    orientations = signs, fits = diagnostics, pivot = pivot,
    centre = centre, scale = spread
  )
  result
}
