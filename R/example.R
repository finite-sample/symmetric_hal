# Rscript R/example.R
source("R/basis.R")
source("R/fit.R")
x <- cbind(1:21, 1:21)
y <- as.numeric(x[, 1] >= 11)
newx <- rbind(c(1, 21), c(21, 1))
signs <- as.matrix(expand.grid(c(1, -1), c(1, -1)))
foldid <- rep(1:7, 3)

fit_example <- function(x, family, tuning) {
  if (tuning == "CV") {
    cv_hal(x, y, c(0.1, 0.05, 0.01), foldid, family = family)
  } else {
    fit_hal(x, y, 0.01, family = family)
  }
}

rows <- list()
for (family in c("gaussian", "binomial")) for (tuning in c("Fixed", "CV")) {
  original <- fit_example(x, family, tuning)
  for (i in seq_len(nrow(signs))) {
    reflected_x <- sweep(x, 2, signs[i, ], `*`)
    fit <- fit_example(reflected_x, family, tuning)
    prediction <- fit$predict(sweep(newx, 2, signs[i, ], `*`))
    residuals <- fit$diagnostics
    if (tuning == "CV") {
      residuals <- c(residuals, unlist(fit$cv_diagnostics[, -(1:2)]))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      family = family, tuning = tuning,
      sign1 = signs[i, 1], sign2 = signs[i, 2],
      lambda = fit$lambda, first_prediction = prediction[1],
      second_prediction = prediction[2],
      test_movement = max(abs(prediction - original$predict(newx))),
      training_movement = max(abs(fit$predict(reflected_x) -
                                    original$predict(x))),
      optimality_residual = max(residuals)
    )
  }
}
results <- do.call(rbind, rows)
stopifnot(max(results$test_movement) < 1e-6,
          max(results$optimality_residual) < 1e-6)
dir.create("data", showWarnings = FALSE)
dir.create("paper/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "data/example.csv", row.names = FALSE)
summary <- aggregate(
  cbind(test_movement, training_movement, optimality_residual) ~
    family + tuning,
  results, max
)
display <- summary
for (column in names(display)[-(1:2)]) {
  display[[column]] <- formatC(display[[column]], format = "e", digits = 2)
}
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Maximum movement across four sign codings in the diagonal",
  "example. Fixed fits use $\\lambda=0.01$; CV uses seven shared folds",
  "and candidates $0.01,0.05,0.1$. Residuals include both optimization",
  "stages and, for CV, all fold fits.}", "\\label{tab:example}",
  knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
               col.names = c("Family", "Tuning", "Test movement",
                             "Training movement", "Residual")),
  "\\end{table}"
), "paper/tables/example.tex")
writeLines(trimws(capture.output(utils::sessionInfo()), which = "right"),
           "data/example_session.txt")
print(summary, row.names = FALSE)
