# Gaussian and binomial HAL with minimum-norm selection of lasso solutions.

hal_inputs <- function(x, y, family, lambda) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) < 1 ||
        nrow(x) < 2 || any(!is.finite(x))) {
    stop("x must be a finite numeric matrix with at least two rows")
  }
  if (!is.numeric(y) || length(y) != nrow(x) || any(!is.finite(y))) {
    stop("y must be finite and have one value per row of x")
  }
  if (!family %in% c("gaussian", "binomial")) stop("unsupported family")
  if (family == "binomial" &&
        (!all(y %in% c(0, 1)) || length(unique(y)) != 2)) {
    stop("binomial training outcomes must contain both 0 and 1")
  }
  if (!is.numeric(lambda) || !length(lambda) ||
        any(!is.finite(lambda) | lambda <= 0)) {
    stop("lambda must contain finite positive values")
  }
}

hal_loss <- function(eta, y, family) {
  if (family == "gaussian") {
    return((eta - y)^2 / 2)
  }
  pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta
}

hal_kkt <- function(design, y, eta, beta, lambda, family) {
  score <- if (family == "gaussian") eta - y else plogis(eta) - y
  gradient <- as.vector(crossprod(design, score)) / length(y)
  residual <- ifelse(beta != 0, abs(gradient + lambda * sign(beta)),
    pmax(abs(gradient) - lambda, 0)
  )
  list(residual = max(abs(mean(score)), residual), gradient = gradient,
       score = score)
}

# Lasso KKT conditions fix zero coefficients outside the equicorrelation set
# and signs inside it. On that face, fixing eta also fixes the L1 norm.
hal_min_norm <- function(design, eta, beta, gradient, lambda,
                         face_tol, score) {
  selected <- which(abs(gradient) >= lambda - min(face_tol, lambda / 2) |
                      beta != 0)
  result <- rep(0, ncol(design))
  if (!length(selected)) {
    return(list(beta = result, kkt = 0, global_kkt = 0))
  }
  signs <- -sign(gradient[selected])
  if (any(signs * beta[selected] < 0)) stop("inconsistent lasso face signs")
  centered <- sweep(
    design[, selected, drop = FALSE], 2,
    colMeans(design[, selected, drop = FALSE]), `-`
  )
  signed <- sweep(centered, 2, signs, `*`)
  decomposition <- svd(signed)
  rank <- sum(decomposition$d > max(dim(signed)) * .Machine$double.eps *
                max(decomposition$d))
  if (!rank) {
    return(list(beta = result, kkt = 0, global_kkt = 0))
  }
  indices <- seq_len(rank)
  equality <- decomposition$v[, indices, drop = FALSE]
  target <- as.vector(crossprod(
    decomposition$u[, indices, drop = FALSE],
    eta - mean(eta)
  )) / decomposition$d[indices]
  p <- length(selected)
  constraints <- cbind(equality, diag(p))
  qp <- quadprog::solve.QP(diag(p), rep(0, p), constraints,
    c(target, rep(0, p)),
    meq = rank
  )
  z <- qp$solution
  # Active nonnegativity constraints fix these coordinates to zero. Preserve
  # that status when roundoff leaves a tiny positive solver value.
  z[qp$iact[qp$iact > rank] - rank] <- 0
  inequality_dual <- qp$Lagrangian[rank + seq_len(p)]
  # quadprog reports unsigned equality multipliers; recover their signs in
  # the orthonormal equality basis before checking stationarity.
  equality_dual <- as.vector(crossprod(equality, z - inequality_dual))
  kkt <- max(
    abs(z - as.vector(equality %*% equality_dual) - inequality_dual),
    abs(as.vector(crossprod(equality, z)) - target),
    pmax(-z, -inequality_dual, 0), abs(z * inequality_dual)
  )
  result[selected] <- signs * pmax(z, 0)
  centered_all <- sweep(design, 2, colMeans(design), `-`)
  equality_vector <- as.vector(decomposition$u[, indices, drop = FALSE] %*%
                                 (equality_dual / decomposition$d[indices]))
  direction <- -(score - mean(score)) / (nrow(design) * lambda)
  direction_error <- 1 - as.vector(crossprod(signed, direction))
  correction <- as.vector(crossprod(equality, direction_error))
  direction <- direction + as.vector(
    decomposition$u[, indices, drop = FALSE] %*%
      (correction / decomposition$d[indices])
  )
  correlation <- as.vector(crossprod(centered_all, direction))
  linear <- -as.vector(crossprod(centered_all, equality_vector))
  inactive <- result == 0
  denominator <- ifelse(linear >= 0, 1 + correlation, 1 - correlation)
  if (any(denominator[inactive] <= 0 & abs(linear[inactive]) > 1e-10)) {
    stop("cannot certify the complete minimum-norm face")
  }
  bound <- ifelse(denominator[inactive] > 0,
                  abs(linear[inactive]) / denominator[inactive], 0)
  multiplier <- max(0, bound)
  secondary_gradient <- linear - multiplier * correlation
  global_kkt <- max(
    ifelse(inactive, pmax(abs(secondary_gradient) - multiplier, 0),
           abs(result + secondary_gradient + multiplier * sign(result))),
    abs(multiplier * (sum(abs(result)) - sum(abs(beta))))
  )
  list(beta = result, kkt = kkt, global_kkt = global_kkt)
}

fit_hal <- function(x, y, lambda, num_knots = 3, symmetric = TRUE,
                    family = "gaussian", solver_tol = 1e-20,
                    certificate_tol = 1e-6, face_tol = 1e-8, knots = NULL,
                    knot_rule = "quantile") {
  hal_inputs(x, y, family, lambda)
  knot_rule <- match.arg(knot_rule, c("quantile", "observed"))
  if (!is.null(knots) && knot_rule != "quantile") {
    stop("supplied knots cannot be combined with an observed knot rule")
  }
  if (length(lambda) != 1L) stop("fixed-lambda fit requires one lambda")
  if (length(num_knots) != 1L || !is.finite(num_knots) ||
        num_knots < 1 || num_knots != as.integer(num_knots)) {
    stop("num_knots must be a positive integer")
  }
  if (length(solver_tol) != 1L || !is.finite(solver_tol) || solver_tol <= 0 ||
        length(certificate_tol) != 1L || !is.finite(certificate_tol) ||
        certificate_tol <= 0 || length(face_tol) != 1L ||
        !is.finite(face_tol) || face_tol <= 0) {
    stop("tolerances must be finite and positive")
  }
  if (is.null(knots)) {
    features <- make_dictionary(x, num_knots, symmetric, knot_rule)
  } else {
    valid <- is.list(knots) && length(knots) == ncol(x) &&
      all(vapply(knots, function(k) {
        is.numeric(k) && length(k) > 0 && all(is.finite(k))
      }, logical(1)))
    if (!valid) stop("knots must contain one finite numeric vector per column")
    features <- dictionary_features(lapply(knots, function(k) sort(unique(k))),
                                    symmetric)
  }
  design <- as.matrix(design_matrix(x, features))
  p <- ncol(design)
  n <- nrow(design)
  varying <- apply(design, 2, function(z) any(z != z[1]))
  if (!any(varying) || (family == "gaussian" && all(y == y[1]))) {
    primary_beta <- rep(0, p)
    primary_eta <- rep(if (family == "gaussian") {
      mean(y)
    } else {
      qlogis(mean(y))
    }, n)
  } else {
    solver_design <- if (p == 1) cbind(design, 0) else design
    path <- lambda
    if (family == "binomial") {
      lambda_max <- max(abs(crossprod(design, y - mean(y)))) / n
      if (lambda_max > lambda) {
        path <- exp(seq(log(lambda_max), log(lambda), length.out = 30))
        path[c(1, length(path))] <- c(lambda_max, lambda)
      }
    }
    primary <- glmnet::glmnet(solver_design, y,
      family = family,
      lambda = path, standardize = FALSE,
      control = list(
        thresh = solver_tol,
        maxit = 1000000
      )
    )
    if (primary$jerr != 0 || length(primary$lambda) != length(path) ||
          abs(tail(primary$lambda, 1) / lambda - 1) >
            10 * .Machine$double.eps) {
      stop("primary ", family, " lasso did not converge at lambda ", lambda)
    }
    coefficients <- as.vector(stats::coef(primary)[, length(primary$lambda)])
    primary_beta <- coefficients[1 + seq_len(p)]
    primary_eta <- coefficients[1] + as.vector(design %*% primary_beta)
  }
  primary_check <- hal_kkt(
    design, y, primary_eta, primary_beta,
    lambda, family
  )
  secondary <- hal_min_norm(
    design, primary_eta, primary_beta,
    primary_check$gradient, lambda,
    max(face_tol, 10 * primary_check$residual), primary_check$score
  )
  coef <- secondary$beta
  b0 <- mean(primary_eta - as.vector(design %*% coef))
  final_eta <- b0 + as.vector(design %*% coef)
  final_check <- hal_kkt(design, y, final_eta, coef, lambda, family)
  budget <- sum(abs(primary_beta))
  primary_objective <- mean(hal_loss(primary_eta, y, family)) +
    lambda * budget
  objective <- mean(hal_loss(final_eta, y, family)) +
    lambda * sum(abs(coef))
  diagnostics <- c(
    primary_kkt = primary_check$residual,
    final_kkt = final_check$residual,
    secondary_kkt = secondary$kkt,
    global_secondary_kkt = secondary$global_kkt,
    prediction_residual = max(abs(final_eta - primary_eta)),
    budget_change = abs(sum(abs(coef)) - budget),
    objective_change = abs(objective - primary_objective)
  )
  if (any(!is.finite(diagnostics)) || max(diagnostics) > certificate_tol) {
    stop("HAL optimality check failed: residual ", max(diagnostics))
  }
  list(
    beta = coef, intercept = b0, lambda = lambda, features = features,
    fitted = if (family == "gaussian") final_eta else plogis(final_eta),
    objective = objective, primary_objective = primary_objective,
    diagnostics = diagnostics, family = family,
    knot_rule = if (is.null(knots)) knot_rule else "supplied",
    predict = function(newx, type = "response") {
      if (!is.matrix(newx) || !is.numeric(newx) ||
            ncol(newx) != ncol(x) || any(!is.finite(newx))) {
        stop("newx must be a finite matrix with the training columns")
      }
      if (!type %in% c("response", "link")) stop("unsupported type")
      linear <- b0 + as.vector(design_matrix(newx, features) %*% coef)
      if (family == "binomial" && type == "response") {
        plogis(linear)
      } else {
        linear
      }
    }
  )
}

cv_hal <- function(x, y, lambda, foldid, num_knots = 3,
                   symmetric = TRUE, family = "gaussian",
                   tie_tol = 1e-8, solver_tol = 1e-20,
                   certificate_tol = 1e-6, knots = NULL,
                   knot_rule = "quantile") {
  hal_inputs(x, y, family, lambda)
  knot_rule <- match.arg(knot_rule, c("quantile", "observed"))
  if (length(foldid) != nrow(x) || !is.numeric(foldid) ||
        any(!is.finite(foldid)) || any(foldid != as.integer(foldid)) ||
        length(unique(foldid)) < 2) {
    stop("supply at least two explicit CV folds")
  }
  if (length(tie_tol) != 1L || !is.finite(tie_tol) || tie_tol < 0) {
    stop("tie_tol must be finite and nonnegative")
  }
  lambda <- sort(unique(lambda), decreasing = TRUE)
  losses <- matrix(NA_real_, nrow(x), length(lambda))
  validation_predictions <- losses
  start <- proc.time()[["elapsed"]]
  certificates <- list()
  for (fold in sort(unique(foldid))) {
    validation <- foldid == fold
    for (j in seq_along(lambda)) {
      fit <- fit_hal(
        x[!validation, , drop = FALSE], y[!validation],
        lambda[j], num_knots, symmetric, family,
        solver_tol, certificate_tol, knots = knots, knot_rule = knot_rule
      )
      eta <- fit$predict(x[validation, , drop = FALSE], type = "link")
      validation_predictions[validation, j] <- eta
      losses[validation, j] <- hal_loss(eta, y[validation], family)
      certificates[[length(certificates) + 1L]] <- data.frame(
        fold = fold, lambda = lambda[j], t(fit$diagnostics), check.names = FALSE
      )
    }
  }
  validation_seconds <- proc.time()[["elapsed"]] - start
  mean_loss <- colMeans(losses)
  best <- which(mean_loss <= min(mean_loss) + tie_tol)[1]
  start <- proc.time()[["elapsed"]]
  fit <- fit_hal(
    x, y, lambda[best], num_knots, symmetric, family,
    solver_tol, certificate_tol, knots = knots, knot_rule = knot_rule
  )
  fit$seconds <- c(validation = validation_seconds,
                   refit = proc.time()[["elapsed"]] - start)
  fit$cv <- data.frame(lambda = lambda, loss = mean_loss)
  fit$validation_predictions <- validation_predictions
  fit$cv_diagnostics <- do.call(rbind, certificates)
  fit$foldid <- foldid
  fit
}

# Each fitter returns list(predict = function(newx), beta, ...). `beta` is
# over the dictionary's features (NULL for tree and forest methods).

fit_cart <- function(x, y) {
  tree <- rpart::rpart(y ~ ., data = data.frame(y = y, x), method = "anova",
                       control = rpart::rpart.control(cp = 0, xval = 10))
  cp <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
  tree <- rpart::prune(tree, cp = cp)
  list(predict = function(newx) unname(predict(tree, data.frame(newx))),
       tree = tree)
}

fit_rf <- function(x, y) {
  forest <- ranger::ranger(y ~ ., data = data.frame(y = y, x),
                           num.threads = 1)
  list(predict = function(newx) {
    predict(forest, data.frame(newx), num.threads = 1)$predictions
  })
}
