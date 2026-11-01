# Enumerate anchors so every original-coordinate prediction is comparable.
anchor_signs <- function(d) {
  as.matrix(expand.grid(rep(list(c(1, -1)), d)))
}

cv_anchors <- function(x, y, lambda, foldid, num_knots = 3,
                       family = "gaussian", tie_tol = 1e-8, knots = NULL) {
  signs <- anchor_signs(ncol(x))
  oriented_knots <- function(i) {
    if (is.null(knots)) return(NULL)
    lapply(seq_len(ncol(x)), function(j) knots[[j]] * signs[i, j])
  }
  fits <- lapply(seq_len(nrow(signs)), function(i) {
    cv_hal(sweep(x, 2, signs[i, ], `*`), y, lambda, foldid,
           num_knots, symmetric = FALSE, family = family, tie_tol = tie_tol,
           knots = oriented_knots(i))
  })
  lambda <- fits[[1]]$cv$lambda
  validation <- lapply(fits, `[[`, "validation_predictions")
  curves <- do.call(rbind, lapply(seq_along(fits), function(i) {
    data.frame(anchor = i, fits[[i]]$cv)
  }))
  eligible <- curves$loss <= min(curves$loss) + tie_tol
  candidates <- curves[eligible, ]
  selected <- candidates[order(-candidates$lambda, candidates$anchor)[1], ]
  averaged <- Reduce(`+`, lapply(validation, function(eta) {
    if (family == "binomial") plogis(eta) else eta
  })) / length(fits)
  average_loss <- vapply(seq_along(lambda), function(j) {
    if (family == "gaussian") return(mean((averaged[, j] - y)^2) / 2)
    mean(ifelse(y == 1, -log(averaged[, j]), -log1p(-averaged[, j])))
  }, numeric(1))
  if (!is.finite(min(average_loss))) stop("no finite averaged validation loss")
  best_average <- which(average_loss <= min(average_loss) + tie_tol)[1]
  average_lambda <- lambda[best_average]
  refit <- function(anchor, penalty) {
    if (fits[[anchor]]$lambda == penalty) return(fits[[anchor]])
    start <- proc.time()[["elapsed"]]
    fit <- fit_hal(sweep(x, 2, signs[anchor, ], `*`), y, penalty,
                   num_knots, symmetric = FALSE, family = family,
                   knots = oriented_knots(anchor))
    fit$seconds <- c(refit = proc.time()[["elapsed"]] - start)
    fit
  }
  selected_fit <- refit(selected$anchor, selected$lambda)
  average_fits <- lapply(seq_along(fits), refit, penalty = average_lambda)
  validation_seconds <- sum(vapply(fits, function(fit) {
    unname(fit$seconds["validation"])
  }, numeric(1)))
  predict_each <- function(newx, use_fits = fits) {
    predictions <- vapply(seq_along(use_fits), function(i) {
      use_fits[[i]]$predict(sweep(newx, 2, signs[i, ], `*`))
    }, numeric(nrow(newx)))
    matrix(predictions, nrow = nrow(newx), ncol = length(use_fits))
  }
  list(
    signs = signs, fits = fits, cv = curves,
    validation_predictions = validation,
    average_cv = data.frame(lambda = lambda, loss = average_loss),
    average_validation_predictions = averaged,
    selected_anchor = selected$anchor, selected_lambda = selected$lambda,
    average_lambda = average_lambda,
    selected_fit = selected_fit, average_fits = average_fits,
    seconds = c(
      anchored = sum(fits[[1]]$seconds),
      selected = validation_seconds + unname(selected_fit$seconds["refit"]),
      averaged = validation_seconds + sum(vapply(average_fits, function(fit) {
        unname(fit$seconds["refit"])
      }, numeric(1)))
    ),
    predict_each = predict_each,
    predict_selected = function(newx) {
      selected_fit$predict(sweep(newx, 2, signs[selected$anchor, ], `*`))
    },
    predict_average = function(newx) rowMeans(predict_each(newx, average_fits))
  )
}
