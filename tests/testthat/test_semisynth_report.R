test_that("causal simulation summaries retain estimator and pairing", {
  project <- normalizePath(root)
  temporary <- tempfile("semisynth-summary-")
  dir.create(temporary)
  dir.create(file.path(temporary, "data"))
  previous <- setwd(temporary)
  on.exit({
    setwd(previous)
    unlink(temporary, recursive = TRUE)
  })
  estimates <- expand.grid(replicate = 1:3, method = c("oracle", "hal", "shal"),
                           coding = c("original", "continuous_negated"),
                           estimator = c("AIPW raw", "AIPW", "TMLE"),
                           stringsAsFactors = FALSE)
  estimates$att <- 1800 + estimates$replicate - 2 +
    0.5 * (estimates$method == "hal" & estimates$coding == "continuous_negated")
  estimates$truth <- 1800
  estimates$se <- 1
  estimates$bound_se <- 2
  estimates$nuisance_warnings <- as.integer(estimates$replicate == 3)
  estimates$nuisance_seconds <- 2
  estimates$targeting_seconds <- 1
  estimates$q_targeting_floor <- 1e-12
  estimates$seconds <- 3
  estimates$n <- 100
  estimates$prepared_clipped_q0 <- as.integer(estimates$replicate == 2)
  estimates$prepared_clipped_q1 <- 0
  estimates$prepared_clipped_g <- 0
  estimates$valid_targeting <- TRUE
  for (metric in c("q0_risk", "q0_treated_risk", "q1_treated_risk", "g_risk",
                   "raw_drift", "prepared_drift")) estimates[[metric]] <- 0
  write.csv(estimates, "data/semisynth_att.csv", row.names = FALSE)
  capture.output(source(file.path(project, "R", "semisynth_report.R"),
                        local = new.env(parent = globalenv())))
  summary <- read.csv("data/semisynth_summary.csv")
  movement <- read.csv("data/semisynth_movement.csv")
  oracle <- subset(summary, method == "oracle")
  expect_equal(nrow(summary), 18)
  expect_equal(oracle$bias, rep(0, 6))
  expect_equal(oracle$bias_mcse, rep(1 / sqrt(3), 6))
  expect_equal(oracle$rmse, rep(sqrt(2 / 3), 6))
  expect_equal(oracle$rms_se, rep(1, 6))
  expect_equal(oracle$coverage, rep(1, 6))
  expect_true(all(oracle$coverage_lower < 1))
  expect_equal(subset(movement, method == "hal")$mean_movement, rep(0.5, 3))
  expect_equal(subset(movement, method == "oracle")$mean_movement, rep(0, 3))
  comparisons <- read.csv("data/semisynth_comparisons.csv")
  original <- subset(comparisons, coding == "original")
  expect_equal(original$rmse_ratio, rep(1, 3))
  expect_equal(original$log_rmse_ratio_mcse, rep(0, 3))
  expect_equal(original$mse_difference_mcse, rep(0, 3))
  reflected <- subset(comparisons, coding == "continuous_negated")
  expect_equal(reflected$mse_difference, rep(-0.25, 3))
  expect_equal(reflected$mse_difference_mcse, rep(1 / sqrt(3), 3))
  expect_equal(reflected$rmse_ratio, rep(sqrt(8 / 11), 3))
  expected_mcse <- sd(c(27 / 22, -3 / 11, -21 / 22)) / (2 * sqrt(3))
  expect_equal(reflected$log_rmse_ratio_mcse, rep(expected_mcse, 3))
  expect_equal(reflected$rmse_ratio_lower,
               sqrt(8 / 11) * exp(-qt(0.975, 2) * expected_mcse) * rep(1, 3))
  diagnostics <- read.csv("data/semisynth_diagnostics.csv")
  expect_equal(diagnostics$replicates, rep(3, 6))
  expect_equal(diagnostics$prepared_clipped_q0_replicates, rep(1, 6))
  expect_equal(diagnostics$prepared_clipped_q0_mean, rep(1 / 3, 6))
  expect_equal(diagnostics$nuisance_warning_replicates, rep(1, 6))

  report <- function(input) {
    write.csv(input, "data/semisynth_att.csv", row.names = FALSE)
    capture.output(source(file.path(project, "R", "semisynth_report.R"),
                          local = new.env(parent = globalenv())))
  }
  expect_error(report(estimates[-1, ]), "nrow\\(estimates\\)")
  changed <- estimates
  changed$prepared_clipped_q0[1] <- 10
  expect_error(report(changed), "anyDuplicated")
  changed <- estimates
  changed$truth[1] <- 1900
  expect_error(report(changed), "unique\\(estimates\\$truth\\)")
  changed <- estimates
  changed$att <- changed$truth
  report(changed)
  zero <- read.csv("data/semisynth_comparisons.csv")
  expect_true(all(is.na(zero$rmse_ratio)))
  expect_equal(zero$mse_difference, rep(0, 6))
  changed <- estimates
  changed$att[changed$method == "shal"] <- 1800
  report(changed)
  zero <- read.csv("data/semisynth_comparisons.csv")
  expect_equal(zero$rmse_ratio, rep(0, 6))
  expect_true(all(is.na(zero$log_rmse_ratio_mcse)))
  expect_true(all(is.na(zero$rmse_ratio_lower)))
  expect_true(all(is.na(zero$rmse_ratio_upper)))
})
