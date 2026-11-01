# Rscript R/causal_control_report.R [completed_nuisance_dir ...]
for (file in c("att", "causal", "runs")) {
  source(file.path("R", paste0(file, ".R")))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  runs <- unlist(lapply(args, read_run), recursive = FALSE)
  jobs <- do.call(rbind, lapply(runs, function(run) {
    data.frame(method = run$method, coding = run$coding)
  }))
  stopifnot(!anyDuplicated(jobs))
  reference <- runs[[1]]$predictions
  estimates <- lapply(runs, function(run) {
    p <- run$predictions
    stopifnot(all(vapply(c("row", "fold", "a", "y"), function(column) {
      identical(p[[column]], reference[[column]])
    }, logical(1))))
    comparison <- compare_att(p, steps = c(1e-5, 1e-6),
                              outcome_bounds = range(reference$y))
    data.frame(method = run$method, coding = run$coding,
               nuisance_warnings = length(run$warnings), comparison)
  })
  estimates <- do.call(rbind, estimates)
  final <- subset(estimates, estimator != "TMLE" | step == 1e-6)
  stopifnot(all(final$valid_targeting))
  diagnostics <- do.call(rbind, lapply(runs, function(run) {
    do.call(rbind, lapply(run$diagnostics, function(d) {
      fits <- if (is.null(d$fit$fits)) list(d$fit) else d$fit$fits
      do.call(rbind, lapply(seq_along(fits), function(i) {
        fit <- fits[[i]]
        grid <- fit$lambda_grid
        if (is.null(grid) && is.data.frame(fit$cv)) grid <- fit$cv$lambda
        data.frame(
          method = run$method, coding = run$coding, fold = d$fold,
          component = d$component, n = d$n, seconds = d$seconds,
          orientation = i, lambda = fit$lambda, basis_size = fit$basis_size,
          path_length = length(grid), selected_index = match(fit$lambda, grid),
          selected_at_minimum = if (length(grid)) {
            fit$lambda == min(grid)
          } else {
            NA
          },
          component_warnings = length(d$warnings),
          warning_messages = paste(unique(d$warnings), collapse = " | ")
        )
      }))
    }))
  }))
} else {
  estimates <- read.csv("data/causal_control_estimates.csv")
  final <- subset(estimates, estimator != "TMLE" | step == 1e-6)
  stopifnot(all(final$valid_targeting))
}
expected_floor <- formals(prepare_att)$q_targeting_floor
if (!"q_targeting_floor" %in% names(estimates) ||
      anyNA(estimates$q_targeting_floor) ||
      any(estimates$q_targeting_floor != expected_floor)) {
  stop("incompatible q_targeting_floor; recompute from saved predictions")
}
control_methods <- c("hal9001_tight", "hal9001_scaled", "hal9001_scaled_tight")
controls <- subset(final, method %in% control_methods)
canonical <- subset(final, method == "hal9001_canonical")
require_grid <- function(data, methods, label) {
  expected <- expand.grid(
    method = methods,
    coding = c("original", "continuous_negated", "units"),
    estimator = c("AIPW raw", "AIPW", "TMLE"),
    stringsAsFactors = FALSE
  )
  key <- function(x) paste(x$method, x$coding, x$estimator, sep = "/")
  if (anyDuplicated(key(data)) ||
        !setequal(key(data), key(expected))) {
    stop(label, " table requires all methods, codings, and estimators")
  }
}
require_grid(controls, control_methods, "control")
require_grid(canonical, "hal9001_canonical", "canonical")
scientific <- function(x) {
  if (x == 0) return("0")
  parts <- strsplit(sprintf("%.2e", x), "e", fixed = TRUE)[[1]]
  sprintf("%s\\times10^{%d}", parts[1], as.integer(parts[2]))
}
if (length(args)) {
  readr::write_csv(estimates, "data/causal_control_estimates.csv")
  readr::write_csv(diagnostics, "data/causal_control_fits.csv")
}
if (nrow(controls)) {
  table <- controls[controls$estimator == "TMLE", ]
  codings <- c("original", "continuous_negated", "units")
  table <- table[order(match(table$method, control_methods),
                       match(table$coding, codings)), ]
  display <- table[c("method", "coding", "att", "se", "nuisance_warnings")]
  display$method <- c("Tighter solver", "Standardized", "Both controls")[
    match(table$method, control_methods)
  ]
  display$coding <- c("Original", "Signs reversed", "Different units")[
    match(table$coding, codings)
  ]
  writeLines(c(
    "\\begin{table}[ht]", "\\centering", "\\small",
    "\\caption{Solver and scaling controls for the first-order nuisance",
    "learner. Estimates and empirical influence-function SEs are in dollars.",
    "Tighter fits use convergence threshold $10^{-12}$ and iteration limit",
    "$10^6$; standardization uses each outer training sample. All codings",
    "share observations, folds, and targeting settings. Every reported",
    "TMLE passes the score and step-refinement checks. Warnings count",
    "nuisance-fit messages: passing the targeting checks does not establish",
    "that every nuisance validation path completed.}",
    "\\label{tab:causal-controls}",
    knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
                 digits = c(0, 0, 0, 0, 0),
                 col.names = c("Control", "Coding", "Estimate", "SE",
                               "Warnings"),
                 linesep = c("", "", "\\addlinespace"), escape = FALSE),
    "\\end{table}"
  ), "paper/tables/causal_controls.tex")
  scaled <- subset(table, method != "hal9001_tight")
  original <- subset(scaled, coding == "original")
  units <- subset(scaled, coding == "units")
  original <- original[match(units$method, original$method), ]
  writeLines(c(
    sprintf("\\newcommand{\\StandardizedUnitATTMovement}{%s}",
            scientific(max(abs(units$att - original$att)))),
    sprintf("\\newcommand{\\StandardizedUnitSEMovement}{%s}",
            scientific(max(abs(units$se - original$se))))
  ), "paper/tables/causal_control_macros.tex")
}
if (nrow(canonical)) {
  codings <- c("original", "continuous_negated", "units")
  stopifnot(setequal(canonical$coding, codings))
  table <- canonical[canonical$estimator == "TMLE", ]
  table <- table[match(codings, table$coding), ]
  display <- table[c("coding", "att", "se", "nuisance_warnings")]
  display$coding <- c("Original", "Signs reversed", "Different units")
  writeLines(c(
    "\\begin{table}[ht]", "\\centering", "\\small",
    "\\caption{ATT TMLE after standardization and canonical orientation",
    "of the first-order nuisance learner. Preprocessing uses each outer",
    "training sample; ambiguous orientations are averaged after separate",
    "package-default tuning. All codings use common folds and targeting",
    "settings. Warnings count nuisance-fit messages and are retained even",
    "when predictions agree. SEs are empirical influence-function SEs.}",
    "\\label{tab:causal-canonical}",
    knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
                 digits = c(0, 0, 0, 0),
                 col.names = c("Coding", "Estimate", "SE", "Warnings"),
                 escape = FALSE), "\\end{table}"
  ), "paper/tables/causal_canonical.tex")
  original <- table$att[1]
  writeLines(c(
    sprintf("\\newcommand{\\CanonicalATT}{%.0f}", original),
    sprintf("\\newcommand{\\CanonicalSignMovement}{%s}",
            scientific(abs(table$att[2] - original))),
    sprintf("\\newcommand{\\CanonicalUnitMovement}{%s}",
            scientific(abs(table$att[3] - original))),
    sprintf("\\newcommand{\\CanonicalWarnings}{%d}",
            table$nuisance_warnings[1])
  ), "paper/tables/causal_canonical_macros.tex")
}
print(final[c("method", "coding", "estimator", "att", "se",
              "nuisance_warnings")], row.names = FALSE)
