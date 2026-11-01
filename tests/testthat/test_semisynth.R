test_that("semi-synthetic truth and population efficiency benchmark agree", {
  source(file.path(root, "R", "semisynth.R"), local = TRUE)
  dgp <- semisynth_dgp(read.csv(file.path(root, "data", "lalonde_dw.csv")))
  population <- dgp$population
  expect_equal(population$true_q1 - population$true_q0,
               rep(1800, nrow(population)))
  expect_true(all(population$true_g >= 0 & population$true_g <= 0.95))
  g <- population$true_g
  p <- mean(g)
  conditional_variance <- g * (dgp$sigma / p)^2 +
    (1 - g) * (dgp$sigma * g / ((1 - g) * p))^2
  expect_equal(dgp$efficiency_variance, mean(conditional_variance))
  first <- semisynth_draw(dgp, 1, 400)
  second <- semisynth_draw(dgp, 2, 400)
  expect_identical(first, semisynth_draw(dgp, 1, 400))
  expect_equal(first$bound_se, second$bound_se)
  expect_false(identical(first$truth$source_row, second$truth$source_row))
  expect_equal(first$truth$true_q0,
               population$true_q0[first$truth$population_row])
})

test_that("oracle initial nuisances have zero raw risk and drift", {
  source(file.path(root, "R", "semisynth.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  dgp <- semisynth_dgp(read.csv(file.path(root, "data", "lalonde_dw.csv")))
  sim <- semisynth_draw(dgp, 1, 400)
  p <- data.frame(a = sim$a, y = sim$y, q0 = sim$truth$true_q0,
                  q1 = sim$truth$true_q1, g = sim$truth$true_g)
  loss <- semisynth_nuisance_loss(p, sim$truth)
  expect_equal(as.numeric(loss[c("q0_risk", "q0_treated_risk",
                                 "q1_treated_risk", "g_risk", "raw_drift")]),
               rep(0, 5))
})

test_that("nuisance drift agrees with the fixed-function ATT score mean", {
  source(file.path(root, "R", "semisynth.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  truth <- data.frame(true_q0 = c(2, 4, 1, 3), true_g = c(0.2, 0.4, 0.6, 0.8))
  effect <- 2
  truth$true_q1 <- truth$true_q0 + effect
  p <- data.frame(a = c(0, 1, 0, 1), y = c(1, 5, 2, 7),
                  q0 = c(0, 5, -2, 6), q1 = c(2, 7, 0, 8),
                  g = c(0.1, 0.7, 0.9, 0.999))
  loss <- semisynth_nuisance_loss(p, truth)
  score_mean <- function(q, g) {
    treated <- truth$true_g * (truth$true_q1 - q - effect)
    control <- (1 - truth$true_g) * g / (1 - g) * (truth$true_q0 - q)
    mean(treated - control) / mean(truth$true_g)
  }
  expect_equal(loss$q0_risk, 5.75)
  expect_equal(loss$q0_treated_risk, 6.9)
  expect_equal(loss$q1_treated_risk, 6.9)
  expect_equal(loss$raw_drift, 443 / 18)
  expect_equal(loss$raw_drift, score_mean(p$q0, pmin(p$g, 0.99)))
  prepared <- prepare_att(p)
  q <- prepared$outcome_bounds[1] + prepared$spread * prepared$q[, "Q0W"]
  expect_equal(loss$prepared_drift, score_mean(q, prepared$g))
  expect_gt(abs(loss$prepared_drift - loss$raw_drift), 1)
})
