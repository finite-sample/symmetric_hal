# ATT targeting with a fixed sample and externally cross-fitted nuisances.
prepare_att <- function(predictions, g_bounds = c(0.01, 0.99),
                        q_bound = 1e-6, outcome_bounds = range(predictions$y),
                        targeting_floor = 1e-8, q_targeting_floor = 1e-12) {
  p <- predictions
  stopifnot(is.finite(targeting_floor), targeting_floor > 0,
            targeting_floor < 0.5)
  stopifnot(all(is.finite(as.matrix(p[c("a", "y", "q0", "q1", "g")]))),
            all(p$a %in% c(0, 1)), all(c(0, 1) %in% p$a),
            length(g_bounds) == 2, 0 < g_bounds[1], g_bounds[2] < 1,
            g_bounds[1] < g_bounds[2], 0 < q_bound, q_bound < 0.5)
  stopifnot(length(q_targeting_floor) == 1, is.finite(q_targeting_floor),
            q_targeting_floor > 0, q_targeting_floor < q_bound,
            1 - q_targeting_floor < 1)
  spread <- diff(outcome_bounds)
  if (length(outcome_bounds) != 2 || !is.finite(spread) || spread <= 0 ||
        any(p$y < outcome_bounds[1] | p$y > outcome_bounds[2])) {
    stop("outcome bounds must contain y and have positive finite width")
  }
  y <- (p$y - outcome_bounds[1]) / spread
  q0 <- (p$q0 - outcome_bounds[1]) / spread
  q1 <- (p$q1 - outcome_bounds[1]) / spread
  clamp_q <- function(q) pmin(pmax(q, q_bound), 1 - q_bound)
  q <- cbind(QAW = (1 - p$a) * clamp_q(q0) + p$a * clamp_q(q1),
             Q0W = clamp_q(q0), Q1W = clamp_q(q1))
  clipped_g <- pmin(pmax(p$g, g_bounds[1]), g_bounds[2])
  calibration <- stats::glm(p$a ~ 1, offset = qlogis(clipped_g),
                            family = stats::binomial())
  g <- as.vector(stats::fitted(calibration))
  if (abs(mean(g) - mean(p$a)) > 1e-8) stop("propensity calibration failed")
  list(a = p$a, y = y, q = q, g = g, spread = spread,
       outcome_bounds = outcome_bounds, g_bounds = g_bounds, q_bound = q_bound,
       targeting_floor = targeting_floor, q_targeting_floor = q_targeting_floor,
       calibration_intercept = unname(stats::coef(calibration)[1]),
       mean_g_before = mean(clipped_g), mean_g_after = mean(g),
       clipped_q0 = sum(q0 != q[, "Q0W"]),
       clipped_q1 = sum(q1 != q[, "Q1W"]), clipped_g = sum(clipped_g != p$g))
}

compare_aipw_codings <- function(original, recoded,
                                 outcome_bounds = range(original$y)) {
  required <- c("row", "fold", "a", "y", "q0", "q1", "g")
  valid_rows <- function(p) {
    all(required %in% names(p)) && nrow(p) >= 2 &&
      !anyNA(p$row) && !anyDuplicated(p$row)
  }
  if (!valid_rows(original) || !valid_rows(recoded)) {
    stop("coding comparisons require unique rows and all nuisance columns")
  }
  original <- original[order(original$row), ]
  recoded <- recoded[order(recoded$row), ]
  same <- vapply(c("row", "fold", "a", "y"), function(column) {
    identical(original[[column]], recoded[[column]])
  }, logical(1))
  if (!all(same)) {
    stop("codings must share rows, folds, treatments and outcomes")
  }
  effective <- function(p) {
    prepared <- prepare_att(p, outcome_bounds = outcome_bounds)
    list(
      "AIPW raw" = list(q0 = p$q0, g = pmin(p$g, 0.99)),
      "AIPW" = list(
        q0 = outcome_bounds[1] + prepared$spread * prepared$q[, "Q0W"],
        g = prepared$g
      )
    )
  }
  first <- effective(original)
  second <- effective(recoded)
  do.call(rbind, lapply(names(first), function(estimator) {
    before <- first[[estimator]]
    after <- second[[estimator]]
    contrast <- aipw_coding_contrast(original$a, original$y,
                                     before$q0, before$g, after$q0, after$g)
    reverse <- aipw_coding_contrast(original$a, original$y,
                                    after$q0, after$g, before$q0, before$g)
    data.frame(
      estimator = estimator, contrast$summary,
      reverse_outcome_component = -reverse$summary$outcome_component,
      reverse_propensity_component = -reverse$summary$propensity_component,
      max_outcome_change = max(abs(after$q0 - before$q0)),
      max_propensity_change = max(abs(after$g - before$g))
    )
  }))
}

target_att <- function(prepared, step = 1e-4, max_iter = 1000000) {
  p <- prepared
  n <- length(p$y)
  fit <- tmle::oneStepATT(
    Y = p$y, A = p$a, Delta = rep(1, n), Q = p$q,
    g1W = p$g, pDelta1 = matrix(1, n, 2),
    depsilon = step, max_iter = max_iter,
    gbounds = c(p$targeting_floor, 1 - p$targeting_floor),
    Qbounds = c(p$q_targeting_floor, 1 - p$q_targeting_floor),
    obsWeights = rep(1, n)
  )
  att <- fit$psi * p$spread
  ic <- fit$IC * p$spread
  if (any(!is.finite(c(att, ic)))) stop("nonfinite TMLE output")
  list(att = att, se = stats::sd(ic) / sqrt(n), ic = ic,
       converged = fit$conv, mean_ic = mean(ic), step = step)
}

compare_att <- function(predictions, steps = c(1e-4, 1e-5, 1e-6),
                        g_bounds = c(0.01, 0.99), q_bound = 1e-6,
                        outcome_bounds = range(predictions$y),
                        q_targeting_floor = 1e-12) {
  stopifnot(length(unique(steps)) >= 2, all(is.finite(steps)), all(steps > 0))
  p <- prepare_att(predictions, g_bounds, q_bound, outcome_bounds,
                   q_targeting_floor = q_targeting_floor)
  raw <- aipw_att(predictions$a, predictions$y, predictions$q0,
                  pmin(predictions$g, g_bounds[2]))
  common <- aipw_att(p$a, p$y, p$q[, "Q0W"], p$g)
  rows <- list(
    data.frame(estimator = "AIPW raw", step = NA_real_, att = raw$att,
               se = raw$se, mean_ic = mean(raw$ic), converged = TRUE),
    data.frame(estimator = "AIPW", step = NA_real_,
               att = common$att * p$spread, se = common$se * p$spread,
               mean_ic = mean(common$ic) * p$spread, converged = TRUE)
  )
  for (step in steps) {
    fit <- target_att(p, step)
    rows[[length(rows) + 1]] <- data.frame(
      estimator = "TMLE", step = step, att = fit$att, se = fit$se,
      mean_ic = fit$mean_ic, converged = fit$converged
    )
  }
  result <- do.call(rbind, rows)
  result$clipped_q0 <- p$clipped_q0
  result$clipped_q1 <- p$clipped_q1
  result$clipped_g <- p$clipped_g
  result$targeting_residual_se <- abs(result$mean_ic) / result$se
  result$calibration_intercept <- p$calibration_intercept
  result$mean_g_before <- p$mean_g_before
  result$mean_g_after <- p$mean_g_after
  result$targeting_floor <- p$targeting_floor
  result$q_targeting_floor <- p$q_targeting_floor
  result$valid_targeting <- result$converged &
    result$targeting_residual_se <= 0.01
  targets <- which(result$estimator == "TMLE")
  result$refinement_att_se <- result$refinement_se_relative <- NA_real_
  if (length(targets) >= 2) {
    final <- targets[which.min(result$step[targets])]
    previous <- targets[order(result$step[targets])[2]]
    result$refinement_att_se[final] <-
      abs(result$att[final] - result$att[previous]) / result$se[final]
    result$refinement_se_relative[final] <-
      abs(result$se[final] - result$se[previous]) / result$se[final]
    result$valid_targeting[final] <- result$valid_targeting[final] &&
      result$refinement_att_se[final] <= 0.01 &&
      result$refinement_se_relative[final] <= 0.01
  }
  result
}
