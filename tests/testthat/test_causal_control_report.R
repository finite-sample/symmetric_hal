test_that("control reports require complete grids and current targeting", {
  project <- normalizePath(root)
  temporary <- tempfile("control-report-")
  dir.create(temporary)
  for (directory in c("R", "data", "paper", "paper/tables")) {
    dir.create(file.path(temporary, directory))
  }
  file.copy(file.path(project, "R", c("att.R", "causal.R", "runs.R")),
            file.path(temporary, "R"))
  previous <- setwd(temporary)
  on.exit({
    setwd(previous)
    unlink(temporary, recursive = TRUE)
  })
  estimates <- expand.grid(
    coding = c("original", "continuous_negated", "units"),
    estimator = c("AIPW raw", "AIPW", "TMLE"), stringsAsFactors = FALSE
  )
  estimates$method <- "hal9001_canonical"
  estimates$step <- ifelse(estimates$estimator == "TMLE", 1e-6, NA_real_)
  estimates$att <- 5
  estimates$se <- 1
  estimates$nuisance_warnings <- 0
  estimates$valid_targeting <- TRUE
  estimates$q_targeting_floor <- 1e-12
  environment <- new.env(parent = globalenv())
  environment$commandArgs <- function(...) character()
  report <- function(input) {
    write.csv(input, "data/causal_control_estimates.csv", row.names = FALSE)
    capture.output(source(file.path(project, "R", "causal_control_report.R"),
                          local = environment))
  }
  expect_error(report(estimates), "control table requires all methods")
  controls <- do.call(rbind, lapply(
    c("hal9001_tight", "hal9001_scaled", "hal9001_scaled_tight"),
    function(method) {
      result <- estimates
      result$method <- method
      result
    }
  ))
  controls$att <- seq_len(nrow(controls)) * 100
  controls$se <- seq_len(nrow(controls)) * 10
  controls$nuisance_warnings <- seq_len(nrow(controls))
  expect_error(report(controls), "canonical table requires all methods")
  expect_error(report(controls[FALSE, ]), "control table requires all methods")
  shuffled <- controls[rev(seq_len(nrow(controls))), ]
  expect_no_error(report(rbind(shuffled, estimates)))
  table <- readLines("paper/tables/causal_controls.tex")
  rows <- grep("^(Tighter solver|Standardized|Both controls)", table,
               value = TRUE)
  expect_length(rows, 9)
  cells <- trimws(strsplit(rows, "&", fixed = TRUE)[[1]])
  expect_equal(cells[1:4], c("Tighter solver", "Original", "700", "70"))
  expect_match(cells[5], "^7")
  cells <- trimws(strsplit(rows, "&", fixed = TRUE)[[9]])
  expect_equal(cells[1:4], c("Both controls", "Different units", "2700", "270"))
  expect_match(cells[5], "^27")
  macros <- readLines("paper/tables/causal_control_macros.tex")
  expect_true(any(grepl(
    "\\newcommand{\\StandardizedUnitATTMovement}{2.00\\times10^{2}}",
    macros, fixed = TRUE
  )))
  expect_true(any(grepl(
    "\\newcommand{\\StandardizedUnitSEMovement}{2.00\\times10^{1}}",
    macros, fixed = TRUE
  )))
  expect_error(report(controls[-1, ]), "requires all methods")
  expect_error(report(rbind(controls, controls[1, ])), "requires all methods")
  expect_error(report(subset(controls, method != "hal9001_scaled_tight")),
               "requires all methods")
  expect_error(report(rbind(controls, estimates[-1, ])),
               "canonical table requires all methods")
  expect_error(report(rbind(controls, estimates, estimates[1, ])),
               "canonical table requires all methods")
  expect_no_error(report(rbind(controls, estimates)))
  original_data <- readLines("data/causal_control_estimates.csv")
  table_paths <- list.files("paper/tables", full.names = TRUE)
  original_tables <- tools::md5sum(table_paths)
  environment$commandArgs <- function(...) "checkpoint"
  environment$read_run <- function(...) {
    lapply(c("original", "continuous_negated", "units"), function(coding) {
      list(method = "hal9001_tight", coding = coding, warnings = character(),
           predictions = data.frame(row = 1:4, fold = c(1, 2, 1, 2),
                                    a = c(0, 1, 0, 1), y = 1:4),
           diagnostics = list(list(fold = 1, component = "q0", n = 4,
                                   seconds = 1, warnings = character(),
                                   fit = list(lambda = 0.1, basis_size = 2,
                                              lambda_grid = 0.1))))
    })
  }
  environment$compare_att <- function(...) {
    estimates[estimates$coding == "original",
              setdiff(names(estimates), c("method", "coding",
                                          "nuisance_warnings"))]
  }
  expect_error(source(file.path(project, "R", "causal_control_report.R"),
                      local = environment), "control table requires")
  current_data <- readLines("data/causal_control_estimates.csv")
  expect_identical(current_data, original_data)
  expect_identical(tools::md5sum(table_paths), original_tables)
  expect_false(file.exists("data/causal_control_fits.csv"))
  environment$commandArgs <- function(...) character()
  estimates$q_targeting_floor <- 1e-6
  expect_error(report(estimates), "incompatible q_targeting_floor")
  estimates$q_targeting_floor <- NULL
  expect_error(report(estimates), "incompatible q_targeting_floor")
})
