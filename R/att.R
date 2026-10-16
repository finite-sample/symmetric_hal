# Cross-fitted augmented IPW estimate of the effect on the treated. `outer`
# assigns rows to folds; `inner[[k]]` holds the CV folds for fold k's
# training rows (`all`) and its comparison rows (`control`). The outcome
# regression is fit on comparison units only, and the propensity score is
# capped at 0.99 so that the odds g / (1 - g) stay finite.
cross_fit_att <- function(w, a, y, method, outer, inner) {
  n <- length(y)
  q0 <- g <- numeric(n)
  for (k in seq_along(inner)) {
    train <- outer != k
    control <- train & a == 0
    q0[!train] <- fit_nuisance(method, w[control, ], y[control], "gaussian",
                               inner[[k]]$control,
                               list(w[!train, , drop = FALSE]))[[1]]
    g[!train] <- fit_nuisance(method, w[train, ], a[train], "binomial",
                              inner[[k]]$all,
                              list(w[!train, , drop = FALSE]))[[1]]
  }
  g <- pmin(g, 0.99)
  resid <- y - q0
  psi <- sum(a * resid - (1 - a) * g / (1 - g) * resid) / sum(a)
  ic <- (a * (resid - psi) - (1 - a) * g / (1 - g) * resid) / mean(a)
  tibble(att = psi, se = sd(ic) / sqrt(n), max_g_control = max(g[a == 0]))
}

make_folds <- function(a, n_folds = 5, seed = 1) {
  set.seed(seed)
  outer <- sample(rep_len(seq_len(n_folds), length(a)))
  inner <- lapply(seq_len(n_folds), function(k) {
    set.seed(seed * 100 + k)
    list(all = sample(rep_len(1:10, sum(outer != k))),
         control = sample(rep_len(1:10, sum(outer != k & a == 0))))
  })
  list(outer = outer, inner = inner)
}
