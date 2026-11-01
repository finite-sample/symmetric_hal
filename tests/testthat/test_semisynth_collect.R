test_that("semi-synthetic collection combines only distinct paired jobs", {
  project <- normalizePath(root)
  temporary <- tempfile("semisynth-collection-")
  dir.create(temporary)
  for (directory in c("R", "data", "original", "reflected")) {
    dir.create(file.path(temporary, directory))
  }
  file.copy(file.path(project, "R", c("runs.R", "att.R", "causal.R",
                                      "semisynth_collect.R")),
            file.path(temporary, "R"))
  previous <- setwd(temporary)
  on.exit({
    setwd(previous)
    unlink(temporary, recursive = TRUE)
  })
  run <- list(
    job = data.frame(replicate = 1L, method = "oracle", coding = "original"),
    estimates = data.frame(estimator = c("AIPW raw", "TMLE", "TMLE"),
                           step = c(NA, 1e-5, 1e-6), att = c(1700, 1790, 1800),
                           se = 100, valid_targeting = TRUE,
                           clipped_q0 = 1, clipped_q1 = 2, clipped_g = 3),
    predictions = data.frame(row = 1:4, fold = c(1, 2, 1, 2),
                             a = c(0, 1, 0, 1), y = c(10, 20, 30, 40)),
    truth = data.frame(true_g = rep(0.5, 4)), effect = 1800, bound_se = 100,
    nuisance_loss = data.frame(q0_risk = 0), warnings = character(),
    nuisance_seconds = 1, seconds = 2
  )
  for (directory in c("original", "reflected")) {
    run$job$coding <- directory
    saveRDS(run, file.path(directory, "job.rds"))
    configuration <- list(n = 4, codings = list(original = 1, reflected = -1),
                          dgp = list(effect = 1800))
    saveRDS(list(specification = list(jobs = "job.rds",
                                      configuration = configuration,
                                      inputs = tools::md5sum("R/att.R")),
                 session = sessionInfo()),
            file.path(directory, "manifest.rds"))
  }
  environment <- new.env(parent = globalenv())
  environment$commandArgs <- function(...) c("original", "reflected")
  current_estimates <- run$estimates
  environment$compare_att <- function(predictions, steps) {
    expect_equal(steps, c(1e-5, 1e-6))
    expect_equal(predictions$y, c(10, 20, 30, 40))
    current_estimates
  }
  collect <- function() {
    source(file.path(project, "R", "semisynth_collect.R"), local = environment)
  }
  collect()
  result <- read.csv("data/semisynth_att.csv")
  expect_equal(nrow(result), 4)
  expect_equal(result$att, rep(c(1700, 1800), 2))
  expect_equal(result$prepared_clipped_q0, rep(1, 4))
  expect_equal(result$nuisance_warnings, rep(0, 4))
  expect_equal(result$n, rep(4, 4))
  provenance <- paste(readLines("data/semisynth_session.txt"), collapse = "\n")
  expect_match(provenance, unname(tools::md5sum("R/semisynth_collect.R")),
               fixed = TRUE)
  run$estimates$att <- -99999
  run$estimates$valid_targeting <- FALSE
  saveRDS(run, "reflected/job.rds")
  collect()
  refreshed <- read.csv("data/semisynth_att.csv")
  expect_equal(refreshed$att, result$att)
  expect_true(all(refreshed$valid_targeting))
  expect_true(all(refreshed$targeting_seconds >= 0))
  expect_equal(refreshed$seconds,
               refreshed$nuisance_seconds + refreshed$targeting_seconds)
  run$estimates <- NULL
  saveRDS(run, "reflected/job.rds")
  collect()
  expect_equal(read.csv("data/semisynth_att.csv")$att, result$att)
  current_estimates$valid_targeting <- FALSE
  expect_error(collect(), "targeting checks failed.*1.*oracle.*original")
  current_estimates$valid_targeting <- TRUE
  run$predictions$y[1] <- 11
  saveRDS(run, "reflected/job.rds")
  expect_error(collect(), "all\\(vapply")
  environment$commandArgs <- function(...) c("original", "original")
  expect_error(collect(), "anyDuplicated")
  environment$commandArgs <- function(...) c("original", "reflected")
  manifest <- readRDS("reflected/manifest.rds")
  manifest$specification$configuration$dgp$effect <- 1900
  saveRDS(manifest, "reflected/manifest.rds")
  expect_error(collect(), "all\\(vapply")
})
