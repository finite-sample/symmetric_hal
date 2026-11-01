# Rscript R/causal_sensitivity.R
source("R/att.R")
source("R/causal.R")
predictions <- read.csv("data/causal_nuisance.csv")
original <- subset(predictions, coding == "original")
reflected <- subset(predictions, coding == "continuous_negated")
stopifnot(all(vapply(c("row", "fold", "a", "y"), function(column) {
  identical(original[[column]], reflected[[column]])
}, logical(1))))
grid <- expand.grid(q0 = c(FALSE, TRUE), q1 = c(FALSE, TRUE),
                    g = c(FALSE, TRUE))
results <- lapply(seq_len(nrow(grid)), function(i) {
  p <- original
  for (component in names(grid)) {
    if (grid[i, component]) p[[component]] <- reflected[[component]]
  }
  fit <- compare_att(p, steps = c(1e-5, 1e-6))
  data.frame(study = "components", q0_reflected = grid$q0[i],
             q1_reflected = grid$q1[i], g_reflected = grid$g[i], fit)
})
for (negated in c(FALSE, TRUE)) {
  p <- if (negated) reflected else original
  fit <- compare_att(p, steps = c(1e-5, 1e-6), outcome_bounds = c(0, 200000))
  results[[length(results) + 1]] <- data.frame(
    study = "fixed_bounds", q0_reflected = negated, q1_reflected = negated,
    g_reflected = negated, fit
  )
}
result <- do.call(rbind, results)
final <- subset(result, estimator == "TMLE" & step == 1e-6)
stopifnot(all(final$valid_targeting))
readr::write_csv(result, "data/causal_sensitivity.csv")
display <- final[c("study", "q0_reflected", "q1_reflected", "g_reflected",
                   "att", "se")]
print(display, row.names = FALSE)
display$study <- ifelse(display$study == "components", "Components",
                        "Fixed bounds")
for (column in c("q0_reflected", "q1_reflected", "g_reflected")) {
  display[[column]] <- ifelse(display[[column]], "Yes", "No")
}
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{ATT TMLE sensitivity to nuisance components and outcome bounds.",
  "Each component row substitutes the reflected coding's indicated initial",
  "nuisance predictions into the original coding's predictions, then applies",
  "the same calibration and targeting. These are diagnostic combinations,",
  "not an additive decomposition of the full coding difference.",
  "Fixed-bound rows use the original or reflected nuisances with common",
  "outcome bounds of $[0,200000]$ dollars. All rows pass the convergence,",
  "score, and step-refinement checks.}", "\\label{tab:causal-sensitivity}",
  knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
               digits = 0, col.names = c("Study", "$Q_0$ reversed",
                                         "$Q_1$ reversed", "$g$ reversed",
                                         "ATT", "SE"), escape = FALSE),
  "\\end{table}"
), "paper/tables/causal_sensitivity.tex")
g_only <- subset(final, study == "components" & !q0_reflected &
                   !q1_reflected & g_reflected)
bounds <- subset(final, study == "fixed_bounds")
writeLines(c(
  sprintf("\\newcommand{\\TmlePropensitySwap}{%s}",
          format(round(g_only$att), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\TmleBoundOriginal}{%s}",
          format(round(bounds$att[!bounds$g_reflected]), big.mark = ",",
                 trim = TRUE)),
  sprintf("\\newcommand{\\TmleBoundReflected}{%s}",
          format(round(bounds$att[bounds$g_reflected]), big.mark = ",",
                 trim = TRUE))
), "paper/tables/causal_sensitivity_macros.tex")
