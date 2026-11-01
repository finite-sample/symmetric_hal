test_that("pilot collection preserves paired Monte Carlo uncertainty", {
  project <- normalizePath(root)
  temporary <- tempfile("pilot-collection-")
  dir.create(temporary)
  for (directory in c("R", "data", "run", "canonical")) {
    dir.create(file.path(temporary, directory))
  }
  file.copy(file.path(project, "R", "runs.R"), file.path(temporary, "R"))
  previous <- setwd(temporary)
  on.exit({
    setwd(previous)
    unlink(temporary, recursive = TRUE)
  })
  jobs <- paste0("job_", 1:3, ".rds")
  grid <- data.frame(truth = "upper", d = 2, n = 100, snr = 2,
                     distribution = "independent", knot_rule = "quantile",
                     replicate = 1:3, job = 1:3)
  methods <- c("anchored", "selected", "averaged", "symmetric",
               "oracle_anchor", "cart", "rf")
  for (i in 1:3) {
    risk <- ifelse(methods == "symmetric", 2 * i, i)
    result <- data.frame(method = methods, risk = risk,
                         test_risk = risk + ifelse(methods == "symmetric",
                                                   0.25, 0.5),
                         seconds = 0.1, exact_integration = TRUE,
                         max_symmetric_movement = 0, outcome_sd = 1)
    saveRDS(list(cell = grid[i, ], results = result, foldid = 1:3,
                 reflection_movement = c(0, 1e-12 * i),
                 elapsed_seconds = 1),
            file.path("run", jobs[i]))
    other <- result[1, ]
    other$method <- "canonical"
    other$risk <- i / 2
    saveRDS(list(cell = grid[i, ], results = other, foldid = 1:3,
                 reflection_movement = 1e-13 * i, orientations = 2L,
                 elapsed_seconds = 1), file.path("canonical", jobs[i]))
  }
  saveRDS(list(specification = list(jobs = jobs)), "run/manifest.rds")
  saveRDS(list(grid = grid), "run/configuration.rds")
  environment <- new.env(parent = globalenv())
  environment$commandArgs <- function(...) "run"
  capture.output(source(file.path(project, "R", "pilot_collect.R"),
                        local = environment))
  summary <- read.csv("data/pilot_summary.csv")
  contrasts <- read.csv("data/pilot_comparisons.csv")
  selected <- subset(contrasts, comparator == "selected")
  expect_equal(selected$difference, 2)
  expect_equal(selected$difference_mcse, 1 / sqrt(3))
  expect_equal(selected$risk_ratio, 2)
  expect_equal(selected$log_ratio_mcse, 0)
  expect_equal(selected$ratio_lower, 2)
  expect_equal(selected$ratio_upper, 2)
  expect_equal(subset(summary, method == "symmetric")$risk_mcse, 2 / sqrt(3))
  expect_equal(subset(contrasts, comparator == "rf")$difference, 1.75)
  expect_equal(subset(contrasts, comparator == "rf")$log_ratio_mcse,
               abs(2 / 4.25 - 1 / 2.5) / sqrt(3))
  expect_false(subset(contrasts, comparator == "rf")$exact_integration)
  saveRDS(list(specification = list(jobs = jobs)), "canonical/manifest.rds")
  environment$commandArgs <- function(...) c("run", "canonical")
  capture.output(source(file.path(project, "R", "pilot_collect.R"),
                        local = environment))
  contrasts <- read.csv("data/pilot_comparisons.csv")
  canonical <- subset(contrasts, comparator == "canonical")
  expect_equal(canonical$difference, 3)
  expect_equal(canonical$difference_mcse, 1.5 / sqrt(3))
  expect_equal(canonical$risk_ratio, 4)
  diagnostics <- read.csv("data/pilot_diagnostics.csv")
  expect_equal(subset(diagnostics, method == "symmetric")$
                 max_reflection_movement, 1e-12 * (1:3))
  expect_equal(subset(diagnostics, method == "canonical")$
                 max_reflection_movement, 1e-13 * (1:3))
  expect_equal(subset(diagnostics, method == "canonical")$orientations,
               rep(2L, 3))
  other <- readRDS(file.path("canonical", jobs[1]))
  other$foldid <- 3:1
  saveRDS(other, file.path("canonical", jobs[1]))
  expect_error(capture.output(source(file.path(project, "R", "pilot_collect.R"),
                                     local = environment)), "identical")
})
