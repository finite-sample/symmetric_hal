test_that("ATT TMLE preserves affine outcome scaling on exact nuisances", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  a <- rep(c(0, 1), each = 500)
  y <- c(rep(c(0, 1), c(400, 100)), rep(c(0, 1), c(300, 200)))
  p <- data.frame(a = a, y = y, q0 = 0.2, q1 = 0.4, g = 0.5)
  prepared <- prepare_att(p)
  expect_equal(colnames(prepared$q), c("QAW", "Q0W", "Q1W"))
  expect_equal(prepared$q[, "QAW"], ifelse(a == 1, 0.4, 0.2))
  fit <- target_att(prepared)
  expect_equal(fit$att, 0.2, tolerance = 1e-7)
  expect_lt(abs(fit$mean_ic), 1e-7)
  for (column in c("y", "q0", "q1")) p[[column]] <- 2000 + 7000 * p[[column]]
  scaled <- target_att(prepare_att(p))
  expect_equal(scaled$att, 7000 * fit$att, tolerance = 1e-7)
  expect_equal(scaled$se, 7000 * fit$se, tolerance = 1e-7)
  comparison <- compare_att(p)
  expect_equal(comparison$att, rep(1400, nrow(comparison)), tolerance = 1e-5)
})

test_that("treated outcome fits exclude validation observations", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  x <- matrix(1:40, ncol = 1)
  a <- rep(c(0, 1), 20)
  y <- sin(1:40)
  folds <- make_folds(a)
  fit_nuisance <- function(method, x, y, family, foldid, newx) {
    expect_equal(length(foldid), nrow(x))
    expect_length(intersect(x[, 1], newx[[1]][, 1]), 0)
    lapply(newx, function(z) {
      rep(if (family == "binomial") 0.5 else mean(y), nrow(z))
    })
  }
  p <- cross_fit_nuisance(x, a, y, "hal", folds$outer, folds$inner,
                          include_treated = TRUE)
  expect_true(all(is.finite(p$q0)))
  expect_true(all(is.finite(p$q1)))
  expect_equal(p$row, seq_len(40))
  expect_length(attr(p, "diagnostics"), 15)
})

test_that("ATT calibration solves the treatment prevalence score", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  a <- rep(c(0, 1), each = 500)
  y <- c(rep(c(0, 1), c(400, 100)), rep(c(0, 1), c(300, 200)))
  p <- data.frame(a = a, y = y, q0 = 0.2, q1 = 0.4, g = 0.7)
  prepared <- prepare_att(p)
  expect_equal(prepared$mean_g_before, 0.7)
  expect_equal(prepared$mean_g_after, 0.5, tolerance = 1e-8)
  expect_equal(prepared$calibration_intercept, -qlogis(0.7), tolerance = 1e-8)
  comparison <- compare_att(p)
  expect_true(all(comparison$valid_targeting))
  expect_equal(comparison$att, rep(0.2, 5), tolerance = 1e-6)
  expect_error(compare_att(p, steps = 1e-4))
})

test_that("initial outcome clipping does not constrain ATT targeting", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  predictions <- data.frame(a = c(1, 1, 0, 0), y = c(0.4, 0.1, 0.3, 0.5),
                            q0 = 0.2, q1 = c(0, 0.8, 0.4, 0.4), g = 0.5)
  prepared <- prepare_att(predictions, outcome_bounds = c(0, 1))
  expect_equal(unname(prepared$q[1, "QAW"]), prepared$q_bound)
  constrained <- prepared
  constrained$q_targeting_floor <- prepared$q_bound
  stopped <- target_att(constrained, step = 1e-3)
  expect_true(stopped$converged)
  expect_gt(abs(stopped$mean_ic) / stopped$se, 1)
  corrected <- compare_att(predictions, steps = c(1e-3, 1e-4),
                           outcome_bounds = c(0, 1))
  final <- corrected[corrected$estimator == "TMLE" &
                       !is.na(corrected$step) & corrected$step == 1e-4, ]
  expect_true(final$valid_targeting)
  expect_lt(final$targeting_residual_se, 0.01)
  expect_lt(final$refinement_att_se, 0.01)
  expect_gt(abs(final$att - stopped$att), stopped$se)
  expect_equal(final$q_targeting_floor, 1e-12)
  looser <- target_att(prepare_att(predictions, outcome_bounds = c(0, 1),
                                   q_targeting_floor = 1e-14), step = 1e-4)
  expect_equal(looser$att, final$att, tolerance = 1e-10)
  expect_error(prepare_att(predictions, q_targeting_floor = 1e-6))
  expect_error(prepare_att(predictions, q_targeting_floor = 1e-20))
})
