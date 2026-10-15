# Truths on x ~ U(-1, 1)^d. Cut points sit at -0.4, 0 and 0.4, which are
# near (not exactly on) the 9-knot sample-quantile grid, so every HAL
# variant carries the same small approximation error and CART does not.
# All truths need interactions of at most two variables.
box <- function(x, j, l, sides) {
  side <- function(v, s) {
    switch(s, upper = v >= 0.4, lower = v <= -0.4, two = abs(v) <= 0.4)
  }
  as.numeric(side(x[, j], sides[1]) & side(x[, l], sides[2]))
}

many_boxes_spec <- local({
  set.seed(2026)
  data.frame(j = sample(1:5, 8, TRUE), l = sample(6:10, 8, TRUE),
             sj = sample(c("upper", "lower"), 8, TRUE),
             sl = sample(c("upper", "lower"), 8, TRUE),
             coef = sample(c(-1, 1), 8, TRUE))
})

truths <- list(
  tree = function(x) {
    ifelse(x[, 1] >= 0, ifelse(x[, 2] >= 0.4, 2, 0),
           ifelse(x[, 4] >= -0.4, -1, 1))
  },
  additive = function(x) {
    (x[, 1] >= 0) + 2 * (x[, 2] >= -0.4) - 1.5 * (x[, 3] >= 0.4)
  },
  xor = function(x) sign(x[, 1]) * sign(x[, 2]),
  box_upper = function(x) 2 * box(x, 1, 2, c("upper", "upper")),
  box_mixed = function(x) 2 * box(x, 1, 2, c("upper", "lower")),
  box_lower = function(x) 2 * box(x, 1, 2, c("lower", "lower")),
  box_interior = function(x) 2 * box(x, 1, 2, c("two", "two")),
  many_boxes = function(x) {
    s <- many_boxes_spec
    rowSums(vapply(seq_len(nrow(s)), function(i) {
      s$coef[i] * box(x, s$j[i], s$l[i], c(s$sj[i], s$sl[i]))
    }, numeric(nrow(x))))
  },
  friedman1 = function(x) {
    u <- (x + 1) / 2
    10 * sin(pi * u[, 1] * u[, 2]) + 20 * (u[, 3] - 0.5)^2 + 10 * u[, 4] +
      5 * u[, 5]
  }
)

true_vars <- list(
  tree = c(1, 2, 4), additive = 1:3, xor = 1:2, box_upper = 1:2,
  box_mixed = 1:2, box_lower = 1:2, box_interior = 1:2,
  many_boxes = sort(unique(c(many_boxes_spec$j, many_boxes_spec$l))),
  friedman1 = 1:5
)

draw_x <- function(n, d, rho = 0) {
  z <- matrix(rnorm(n * d), n, d)
  if (rho > 0) {
    z <- z %*% chol(rho + (1 - rho) * diag(d))
  }
  x <- 2 * pnorm(z) - 1
  colnames(x) <- paste0("x", seq_len(d))
  x
}

signal_var <- function(truth, d, rho) {
  set.seed(1)
  var(truths[[truth]](draw_x(2e5, d, rho)))
}

methods_all <- c("cart", "hal", "hal_relaxed", "shal", "shal_relaxed", "rf")

fit_method <- function(method, x, y, dicts) {
  switch(method,
    cart = fit_cart(x, y),
    rf = fit_rf(x, y),
    hal = fit_lasso(dicts$a$design, y, dicts$a$features),
    hal_relaxed = fit_lasso(dicts$a$design, y, dicts$a$features,
                            relax = TRUE),
    shal = fit_lasso(dicts$s$design, y, dicts$s$features),
    shal_relaxed = fit_lasso(dicts$s$design, y, dicts$s$features,
                             relax = TRUE)
  )
}

# Variables used and distinct (variable, cut) splits used, comparable
# across trees and dictionary fits.
fit_structure <- function(method, fit, dicts, x_names) {
  if (method == "cart") {
    frame <- fit$tree$frame
    internal <- frame$var != "<leaf>"
    vars <- unique(as.character(frame$var[internal]))
    return(list(vars = match(vars, x_names), n_splits = sum(internal),
                n_coef = sum(!internal)))
  }
  if (is.null(fit$beta)) {
    return(list(vars = NA, n_splits = NA, n_coef = NA))
  }
  f <- if (method %in% c("hal", "hal_relaxed")) dicts$a$features else
    dicts$s$features
  act <- f[fit$beta != 0, ]
  splits <- unique(c(paste(act$j, act$u),
                     paste(act$l, act$v)[!is.na(act$l)]))
  list(vars = sort(unique(c(act$j, act$l[!is.na(act$l)]))),
       n_splits = length(splits), n_coef = nrow(act))
}

run_replicate <- function(truth, snr, rep, n, d, rho, test_x, methods,
                          num_knots = 9) {
  sigma <- sqrt(signal_var(truth, d, rho) / snr)
  set.seed(rep * 7919 + match(truth, names(truths)) * 101 +
             round(log(snr) * 1000) + n + round(rho * 10))
  x <- draw_x(n, d, rho)
  y <- truths[[truth]](x) + rnorm(n, sd = sigma)
  dicts <- list()
  if (any(grepl("^hal", methods))) {
    fa <- make_dictionary(x, num_knots, FALSE)
    dicts$a <- list(features = fa, design = design_matrix(x, fa))
  }
  if (any(grepl("shal", methods))) {
    fs <- make_dictionary(x, num_knots, TRUE)
    dicts$s <- list(features = fs, design = design_matrix(x, fs))
  }
  f_test <- truths[[truth]](test_x)
  null_risk <- mean((f_test - mean(f_test))^2)
  lapply(methods, function(method) {
    start <- proc.time()[["elapsed"]]
    fit <- fit_method(method, x, y, dicts)
    seconds <- proc.time()[["elapsed"]] - start
    pred <- fit$predict(test_x)
    stopifnot(all(is.finite(pred)))
    risk <- mean((pred - f_test)^2)
    st <- fit_structure(method, fit, dicts, colnames(x))
    tv <- true_vars[[truth]]
    sel <- st$vars
    tp <- length(intersect(sel, tv))
    f1 <- if (anyNA(sel)) NA else 2 * tp / (length(sel) + length(tv))
    list(truth = truth, snr = snr, rep = rep, n = n, d = d, rho = rho,
         method = method, risk = risk, rel_risk = risk / null_risk,
         rte = (risk + sigma^2) / sigma^2,
         pve = 1 - (risk + sigma^2) / (null_risk + sigma^2),
         n_splits = st$n_splits, n_coef = st$n_coef, vars = sel, f1 = f1,
         seconds = seconds, pred = pred)
  })
}
