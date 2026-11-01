semisynth_dgp <- function(data, effect = 1800) {
  data$source_row <- seq_len(nrow(data))
  data <- data[data$group %in% c("treated", "psid"), ]
  data$population_row <- seq_len(nrow(data))
  data$u74 <- as.numeric(data$re74 == 0)
  data$u75 <- as.numeric(data$re75 == 0)
  warnings <- character()
  g_model <- withCallingHandlers(
    stats::glm(treat ~ age + I(age^2) + educ + black + hisp + married +
                 nodegree + re74 + re75 + u74 + u75,
               family = stats::binomial(), data = data),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  q_model <- stats::lm(
    re78 ~ splines::ns(age, 3) + educ + black + hisp + married + nodegree +
      splines::ns(re74, 4) + splines::ns(re75, 4) + u74 + u75,
    data = data[data$treat == 0, ]
  )
  data$true_g <- pmin(stats::predict(g_model, data, type = "response"), 0.95)
  data$true_q0 <- as.vector(stats::predict(q_model, data))
  data$true_q1 <- data$true_q0 + effect
  sigma <- stats::sigma(q_model)
  p <- mean(data$true_g)
  stopifnot(g_model$converged, q_model$rank == length(stats::coef(q_model)),
            all(is.finite(as.matrix(data[c("true_g", "true_q0", "true_q1")]))),
            is.finite(sigma), sigma > 0, p > 0, p < 1)
  list(population = data, effect = effect, sigma = sigma, p = p,
       efficiency_variance = sigma^2 * mean(data$true_g / (1 - data$true_g)) /
         p^2, warnings = warnings, g_coefficients = stats::coef(g_model),
       q_coefficients = stats::coef(q_model))
}

semisynth_draw <- function(dgp, replicate, n = nrow(dgp$population)) {
  set.seed(600000 + replicate)
  indices <- sample.int(nrow(dgp$population), n, replace = TRUE)
  data <- dgp$population[indices, ]
  a <- stats::rbinom(n, 1, data$true_g)
  y <- data$true_q0 + dgp$effect * a + stats::rnorm(n, sd = dgp$sigma)
  columns <- c("age", "educ", "black", "hisp", "married", "nodegree",
               "re74", "re75")
  list(w = as.matrix(data[columns]), a = a, y = y,
       truth = data[c("source_row", "population_row", "true_q0", "true_q1",
                      "true_g")],
       effect = dgp$effect, bound_se = sqrt(dgp$efficiency_variance / n))
}

semisynth_nuisance_loss <- function(predictions, truth) {
  prepared <- prepare_att(predictions)
  raw_g <- pmin(predictions$g, 0.99)
  q0 <- prepared$outcome_bounds[1] + prepared$spread * prepared$q[, "Q0W"]
  weights <- truth$true_g / mean(truth$true_g)
  data.frame(
    q0_risk = mean((predictions$q0 - truth$true_q0)^2),
    q0_treated_risk = mean(weights * (predictions$q0 - truth$true_q0)^2),
    q1_treated_risk = mean(weights * (predictions$q1 - truth$true_q1)^2),
    g_risk = mean((predictions$g - truth$true_g)^2),
    raw_drift = mean((raw_g - truth$true_g) / (1 - raw_g) *
                       (predictions$q0 - truth$true_q0)) / mean(truth$true_g),
    prepared_drift = mean((prepared$g - truth$true_g) / (1 - prepared$g) *
                            (q0 - truth$true_q0)) / mean(truth$true_g)
  )
}
