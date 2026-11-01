test_that("AIPW coding contrast equals direct estimate and influence changes", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  a <- c(1, 0, 1, 0, 0, 1)
  y <- c(3, -2, 7, 4, -1, 5)
  q <- c(1, 0, 2, -1, 3, 2)
  g <- c(0.3, 0.2, 0.5, 0.7, 0, 0.6)
  new_q <- q + c(0.5, -1, 1, 2, -0.5, 0)
  new_g <- c(0.4, 0.3, 0.4, 0.99, 0.1, 0.7)
  first <- aipw_att(a, y, q, g)
  second <- aipw_att(a, y, new_q, new_g)
  contrast <- aipw_coding_contrast(a, y, q, g, new_q, new_g)
  result <- contrast$summary
  expect_equal(result$att_change, second$att - first$att)
  expect_equal(result$outcome_component + result$propensity_component,
               result$att_change)
  expect_equal(contrast$contributions$ic_change, second$ic - first$ic)
  expect_lte(abs(result$att_change), result$absolute_bound)
  expect_lte(abs(result$att_change), result$row_absolute_bound)
  expect_lte(result$row_absolute_bound, result$absolute_bound)
  expect_equal(result$positive_row_change + result$negative_row_change,
               result$att_change)
  expect_equal(result$positive_row_change - result$negative_row_change,
               result$row_absolute_bound)
  before_terms <- (a - (1 - a) * g / (1 - g)) * (y - q)
  after_terms <- (a - (1 - a) * new_g / (1 - new_g)) * (y - new_q)
  expect_equal(result$row_absolute_bound,
               sum(abs(after_terms - before_terms)) / sum(a))
  reverse <- aipw_coding_contrast(a, y, new_q, new_g, q, g)$summary
  expect_equal(reverse$row_absolute_bound, result$row_absolute_bound)
  expect_lte(abs(result$se_change), result$se_change_bound)
  expect_equal(result$minimum_control_overlap, 0.01)
  g_only <- aipw_coding_contrast(a, y, q, g, q, new_g)$summary
  q_only <- aipw_coding_contrast(a, y, q, g, new_q, g)$summary
  expect_equal(g_only$outcome_component, 0)
  expect_equal(q_only$propensity_component, 0)
  expect_error(aipw_coding_contrast(a, y, q, g, new_q[-1], new_g),
               "equal-length")
})

test_that("identical AIPW estimates need not have identical standard errors", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  a <- c(1, 1, 0, 0)
  y <- c(1, -1, 1, -1)
  result <- aipw_coding_contrast(a, y, rep(0, 4), rep(0.5, 4),
                                 rep(1, 4), rep(0.5, 4))$summary
  expect_equal(result$original_att, 0)
  expect_equal(result$recoded_att, 0)
  expect_equal(result$original_se, sqrt(4 / 3))
  expect_equal(result$recoded_se, sqrt(8 / 3))
  expect_gt(abs(result$se_change), 0)
})

test_that("coding comparisons align rows and use effective AIPW inputs", {
  source(file.path(root, "R", "att.R"), local = TRUE)
  source(file.path(root, "R", "causal.R"), local = TRUE)
  original <- data.frame(
    row = 1:6, fold = c(1, 2, 1, 2, 1, 2), a = c(1, 0, 1, 0, 0, 1),
    y = c(3, -2, 7, 4, -1, 5), q0 = c(-20, 0, 20, -1, 3, 2),
    q1 = c(3, 0, 5, 2, 3, 4), g = c(0.001, 0.2, 0.5, 0.7, 0, 0.6)
  )
  recoded <- original
  recoded$q0 <- original$q0 + c(0.5, -1, 1, 2, -0.5, 0)
  recoded$g <- c(0.4, 0.3, 0.4, 0.999, 0.1, 0.7)
  result <- compare_aipw_codings(original, recoded)
  shuffled <- recoded[c(5, 1, 3, 2, 6, 4), ]
  expect_equal(compare_aipw_codings(original, shuffled), result)
  prepared_before <- prepare_att(original)
  prepared_after <- prepare_att(recoded)
  old <- aipw_att(original$a, prepared_before$y,
                  prepared_before$q[, "Q0W"], prepared_before$g)
  new <- aipw_att(recoded$a, prepared_after$y,
                  prepared_after$q[, "Q0W"], prepared_after$g)
  calibrated <- result[result$estimator == "AIPW", ]
  expect_equal(calibrated$att_change,
               (new$att - old$att) * prepared_before$spread)
  expect_equal(calibrated$se_change,
               (new$se - old$se) * prepared_before$spread)
  expect_equal(calibrated$minimum_control_overlap,
               min(1 - c(prepared_before$g[original$a == 0],
                         prepared_after$g[recoded$a == 0])))
  raw <- result[result$estimator == "AIPW raw", ]
  old <- aipw_att(original$a, original$y, original$q0, pmin(original$g, 0.99))
  new <- aipw_att(recoded$a, recoded$y, recoded$q0, pmin(recoded$g, 0.99))
  mixed <- aipw_att(original$a, original$y, recoded$q0, pmin(original$g, 0.99))
  expect_equal(raw$att_change, new$att - old$att)
  expect_equal(raw$reverse_outcome_component, mixed$att - old$att)
  expect_equal(raw$reverse_propensity_component, new$att - mixed$att)
  expect_equal(result$reverse_outcome_component +
                 result$reverse_propensity_component, result$att_change)
  expect_gt(abs(raw$outcome_component - raw$reverse_outcome_component), 0)
  changed <- recoded
  changed$y[1] <- changed$y[1] + 1
  expect_error(compare_aipw_codings(original, changed), "must share")
  changed <- recoded
  changed$fold[1] <- 2
  expect_error(compare_aipw_codings(original, changed), "must share")
  expect_error(compare_aipw_codings(original, recoded[-1, ]), "must share")
  expect_error(compare_aipw_codings(original, recoded[c(1, 1:6), ]),
               "unique rows")
})
