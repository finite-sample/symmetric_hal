test_that("paper tables require the full pilot and preserve cell estimates", {
  project <- normalizePath(root)
  temporary <- tempfile("semisynth-paper-")
  dir.create(temporary)
  dir.create(file.path(temporary, "data"))
  dir.create(file.path(temporary, "R"))
  file.copy(file.path(project, "R", "semisynth_report.R"),
            file.path(temporary, "R"))
  previous <- setwd(temporary)
  on.exit({
    setwd(previous)
    unlink(temporary, recursive = TRUE)
  })
  methods <- c("oracle", "hal", "canonical", "shal", "hal9001",
               "hal9001_canonical")
  input <- expand.grid(
    replicate = 1:20, method = methods,
    coding = c("original", "continuous_negated"),
    estimator = c("AIPW raw", "AIPW", "TMLE"), stringsAsFactors = FALSE
  )
  input$truth <- 1800
  input$att <- 1800 + 100 * (input$replicate - 10.5) +
    1000 * (input$method == "hal" & input$coding == "continuous_negated")
  input$se <- 1000
  input$bound_se <- 1000
  input$n <- 2675
  input$q_targeting_floor <- 1e-12
  input$valid_targeting <- TRUE
  input$nuisance_warnings <- as.integer(input$replicate == 1)
  input$nuisance_seconds <- 2
  input$targeting_seconds <- 1
  input$seconds <- 3
  for (metric in c("q0_risk", "q0_treated_risk", "q1_treated_risk", "g_risk",
                   "raw_drift", "prepared_drift", "prepared_clipped_q0",
                   "prepared_clipped_q1", "prepared_clipped_g")) {
    input[[metric]] <- 0
  }
  report <- function(data) {
    write.csv(data, "data/semisynth_att.csv", row.names = FALSE)
    capture.output(source(file.path(project, "R", "semisynth_paper.R"),
                          local = new.env(parent = globalenv())))
  }
  expect_error(report(subset(input, replicate <= 2)), "replicates 1 through 20")
  expect_false(dir.exists("paper/tables"))
  expect_error(report(subset(input, method != "shal")), "all six methods")
  expect_false(dir.exists("paper/tables"))
  expect_error(report(input[-1, ]), "nrow\\(estimates\\)")
  expect_false(dir.exists("paper/tables"))
  report(input)
  expect_setequal(list.files("paper/tables"), paste0("semisynth_", c(
    "aipw_raw", "aipw", "tmle", "movement", "comparisons", "macros"
  ), ".tex"))
  lines <- readLines("paper/tables/semisynth_tmle.tex")
  row <- lines[grepl("Anchored (0th) & Reversed", lines, fixed = TRUE)]
  expect_length(row, 1)
  expect_match(row, "1.00 (0.13)", fixed = TRUE)
  expect_match(row, "1.15 (0.12)", fixed = TRUE)
  expect_match(row, "100 [83, 100]", fixed = TRUE)
  movement <- readLines("paper/tables/semisynth_movement.tex")
  row <- movement[grepl("Anchored (0th) & TMLE", movement, fixed = TRUE)]
  expect_match(row, "1000.00 (0)", fixed = TRUE)
  expect_match(row, "1000.00 & 0", fixed = TRUE)
  macros <- readLines("paper/tables/semisynth_macros.tex")
  expect_true("\\newcommand{\\CausalPilotOracleBias}{0}" %in% macros)
  expect_true("\\newcommand{\\CausalPilotDraws}{20}" %in% macros)

  input$att[input$method == "shal"] <- 1800
  report(input)
  lines <- readLines("paper/tables/semisynth_tmle.tex")
  row <- lines[grepl("Symmetric (0th) & Original", lines, fixed = TRUE)]
  expect_match(row, "0.00 (--)", fixed = TRUE)
  comparisons <- readLines("paper/tables/semisynth_comparisons.tex")
  row <- comparisons[grepl("0th: symmetric / anchored & TMLE", comparisons,
                           fixed = TRUE)]
  expect_match(row, "0.00 [--]", fixed = TRUE)
})
