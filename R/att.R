# Outer-fold predictions for the outcome regressions and propensity score.
cross_fit_nuisance <- function(w, a, y, method, outer, inner,
                               include_treated = FALSE) {
  n <- length(y)
  q0 <- g <- numeric(n)
  q1 <- rep(NA_real_, n)
  diagnostics <- new.env(parent = emptyenv())
  fit_component <- function(component, rows, outcome, family, foldid, heldout) {
    warnings <- character()
    start <- proc.time()[["elapsed"]]
    result <- withCallingHandlers(
      fit_nuisance(method, w[rows, , drop = FALSE], outcome, family,
                   foldid, list(w[heldout, , drop = FALSE])),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
      }
    )
    diagnostics[[paste(k, component, sep = "_")]] <- list(
      fold = k, component = component, n = sum(rows), warnings = warnings,
      seconds = proc.time()[["elapsed"]] - start,
      fit = attr(result, "diagnostics")
    )
    result[[1]]
  }
  for (k in seq_along(inner)) {
    train <- outer != k
    control <- train & a == 0
    q0[!train] <- fit_component("q0", control, y[control], "gaussian",
                                inner[[k]]$control, !train)
    if (include_treated) {
      treated <- train & a == 1
      q1[!train] <- fit_component("q1", treated, y[treated], "gaussian",
                                  inner[[k]]$treated, !train)
    }
    g[!train] <- fit_component("g", train, a[train], "binomial",
                               inner[[k]]$all, !train)
  }
  if (any(!is.finite(q0)) || any(!is.finite(g)) ||
        any(g < 0 | g > 1)) {
    stop("invalid nuisance predictions")
  }
  if (include_treated && any(!is.finite(q1))) {
    stop("invalid treated outcome predictions")
  }
  result <- data.frame(row = seq_len(n), fold = outer, a = a, y = y,
                       q0 = q0, q1 = q1, g = g)
  attr(result, "diagnostics") <- as.list(diagnostics)
  result
}

aipw_att <- function(a, y, q0, g) {
  n <- length(y)
  if (any(!is.finite(c(y, q0, g))) || any(g < 0 | g >= 1) ||
        !all(a %in% c(0, 1)) || !all(c(0, 1) %in% a)) {
    stop("invalid ATT inputs")
  }
  resid <- y - q0
  psi <- sum(a * resid - (1 - a) * g / (1 - g) * resid) / sum(a)
  ic <- (a * (resid - psi) - (1 - a) * g / (1 - g) * resid) / mean(a)
  se <- sd(ic) / sqrt(n)
  if (!is.finite(psi) || !is.finite(se)) stop("nonfinite ATT estimate")
  list(att = psi, se = se, ic = ic)
}

aipw_coding_contrast <- function(a, y, q0, g, new_q0, new_g) {
  n <- length(y)
  inputs <- list(a, q0, g, new_q0, new_g)
  if (n < 2 || any(lengths(inputs) != n)) {
    stop("coding contrasts require equal-length prediction vectors")
  }
  original <- aipw_att(a, y, q0, g)
  recoded <- aipw_att(a, y, new_q0, new_g)
  treated <- sum(a)
  delta_q <- new_q0 - q0
  delta_g <- new_g - g
  outcome <- (-a + (1 - a) * new_g / (1 - new_g)) * delta_q
  propensity <- -(1 - a) * (y - q0) * delta_g /
    ((1 - g) * (1 - new_g))
  change <- sum(outcome + propensity) / treated
  ic_change <- (outcome + propensity - a * change) / mean(a)
  outcome_bound <- sum((a + (1 - a) * new_g / (1 - new_g)) *
                         abs(delta_q)) / treated
  propensity_bound <- sum(abs(propensity)) / treated
  list(
    summary = data.frame(
      original_att = original$att, recoded_att = recoded$att,
      att_change = change, outcome_component = sum(outcome) / treated,
      propensity_component = sum(propensity) / treated,
      outcome_bound = outcome_bound, propensity_bound = propensity_bound,
      absolute_bound = outcome_bound + propensity_bound,
      row_absolute_bound = sum(abs(outcome + propensity)) / treated,
      positive_row_change = sum(pmax(outcome + propensity, 0)) / treated,
      negative_row_change = sum(pmin(outcome + propensity, 0)) / treated,
      original_se = original$se, recoded_se = recoded$se,
      se_change = recoded$se - original$se,
      se_change_bound = stats::sd(ic_change) / sqrt(n),
      minimum_control_overlap = min(1 - c(g[a == 0], new_g[a == 0]))
    ),
    contributions = data.frame(
      row = seq_len(n), outcome = outcome / treated,
      propensity = propensity / treated, ic_change = ic_change
    )
  )
}

cross_fit_att <- function(w, a, y, method, outer, inner) {
  nuisance <- cross_fit_nuisance(w, a, y, method, outer, inner)
  g <- pmin(nuisance$g, 0.99)
  fit <- aipw_att(a, y, nuisance$q0, g)
  tibble::tibble(att = fit$att, se = fit$se, max_g_control = max(g[a == 0]))
}

make_folds <- function(a, n_folds = 5, seed = 1) {
  set.seed(seed)
  outer <- sample(rep_len(seq_len(n_folds), length(a)))
  inner <- lapply(seq_len(n_folds), function(k) {
    set.seed(seed * 100 + k)
    list(all = sample(rep_len(1:10, sum(outer != k))),
         control = sample(rep_len(1:10, sum(outer != k & a == 0))),
         treated = sample(rep_len(1:10, sum(outer != k & a == 1))))
  })
  list(outer = outer, inner = inner)
}
