# Orient each covariate by its centered third moment; average ambiguous signs.
canonical_signs <- function(x, orientation_tol = 1e-12) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 2 || ncol(x) < 1 ||
        any(!is.finite(x))) stop("x must be a finite numeric matrix")
  if (length(orientation_tol) != 1 || !is.finite(orientation_tol) ||
        orientation_tol < 0) stop("orientation_tol must be nonnegative")
  residual <- sweep(x, 2, apply(x, 2, stats::median), `-`)
  centered <- sweep(residual, 2, colMeans(residual), `-`)
  scale <- apply(abs(centered), 2, max)
  scale[scale == 0] <- 1
  score <- colMeans(sweep(centered, 2, scale, `/`)^3)
  choices <- lapply(score, function(value) {
    if (abs(value) <= orientation_tol) c(1, -1) else sign(value)
  })
  signs <- as.matrix(expand.grid(choices))
  attr(signs, "score") <- score
  signs
}

fit_canonical_hal <- function(x, y, lambda, num_knots = 3,
                              family = "gaussian", knots = NULL,
                              orientation_tol = 1e-12) {
  hal_inputs(x, y, family, lambda)
  signs <- canonical_signs(x, orientation_tol)
  fits <- lapply(seq_len(nrow(signs)), function(i) {
    transformed_knots <- if (!is.null(knots)) {
      lapply(seq_len(ncol(x)), function(j) knots[[j]] * signs[i, j])
    }
    fit_hal(sweep(x, 2, signs[i, ], `*`), y, lambda, num_knots,
            symmetric = FALSE, family = family, knots = transformed_knots)
  })
  list(
    lambda = lambda, signs = signs, fits = fits, family = family,
    diagnostics = do.call(rbind, lapply(fits, `[[`, "diagnostics")),
    predict = function(newx, type = "response") {
      if (!type %in% c("response", "link")) stop("unsupported type")
      predictions <- lapply(seq_along(fits), function(i) {
        fits[[i]]$predict(sweep(newx, 2, signs[i, ], `*`))
      })
      average <- Reduce(`+`, predictions) / length(fits)
      if (family == "binomial" && type == "link") qlogis(average) else average
    }
  )
}

cv_canonical_hal <- function(x, y, lambda, foldid, num_knots = 3,
                             family = "gaussian", tie_tol = 1e-8,
                             knots = NULL, orientation_tol = 1e-12) {
  hal_inputs(x, y, family, lambda)
  if (length(foldid) != nrow(x) || !is.numeric(foldid) ||
        any(!is.finite(foldid)) || any(foldid != as.integer(foldid)) ||
        length(unique(foldid)) < 2) {
    stop("supply at least two explicit CV folds")
  }
  if (length(tie_tol) != 1 || !is.finite(tie_tol) || tie_tol < 0) {
    stop("tie_tol must be finite and nonnegative")
  }
  lambda <- sort(unique(lambda), decreasing = TRUE)
  predictions <- losses <- matrix(NA_real_, nrow(x), length(lambda))
  certificates <- list()
  orientations <- list()
  start <- proc.time()[["elapsed"]]
  for (fold in sort(unique(foldid))) {
    validation <- foldid == fold
    for (j in seq_along(lambda)) {
      fit <- fit_canonical_hal(
        x[!validation, , drop = FALSE], y[!validation], lambda[j],
        num_knots, family, knots, orientation_tol
      )
      eta <- fit$predict(x[validation, , drop = FALSE], type = "link")
      predictions[validation, j] <- eta
      losses[validation, j] <- hal_loss(eta, y[validation], family)
      certificates[[length(certificates) + 1L]] <- data.frame(
        fold = fold, lambda = lambda[j], fit$diagnostics, check.names = FALSE
      )
      orientations[[as.character(fold)]] <- fit$signs
    }
  }
  validation_seconds <- proc.time()[["elapsed"]] - start
  mean_loss <- colMeans(losses)
  best <- which(mean_loss <= min(mean_loss) + tie_tol)[1]
  start <- proc.time()[["elapsed"]]
  fit <- fit_canonical_hal(x, y, lambda[best], num_knots, family, knots,
                           orientation_tol)
  fit$seconds <- c(validation = validation_seconds,
                   refit = proc.time()[["elapsed"]] - start)
  fit$cv <- data.frame(lambda = lambda, loss = mean_loss)
  fit$validation_predictions <- predictions
  fit$cv_diagnostics <- do.call(rbind, certificates)
  fit$cv_signs <- orientations
  fit$foldid <- foldid
  fit
}
