# Rscript R/semisynth_oracle_report.R
summary <- read.csv("data/semisynth_oracle_summary.csv")
draws <- read.csv("data/semisynth_oracle.csv")
stopifnot(setequal(summary$interval, c("se", "conditional_se")),
          nrow(summary) == 2, !anyDuplicated(draws$replicate),
          all(summary$replicates == nrow(draws)),
          length(unique(draws$truth)) == 1)
empirical <- subset(summary, interval == "se")
conditional <- subset(summary, interval == "conditional_se")
display <- summary[c("interval", "mean_se", "rms_se", "coverage",
                     "coverage_lower", "coverage_upper")]
labels <- c(se = "Empirical influence function",
            conditional_se = "Known conditional variance")
display$interval <- labels[display$interval]
for (column in c("coverage", "coverage_lower", "coverage_upper")) {
  display[[column]] <- 100 * display[[column]]
}
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Oracle AIPW intervals in the known-effect population.",
  "Both rows use the same estimates with true nuisance functions; only the",
  "standard error changes. Mean and root-mean-square (RMS) standard errors",
  "are in outcome units. Coverage and its exact binomial Monte Carlo interval",
  "are percentages. The known conditional variance is available for this",
  "Gaussian oracle calculation, not for fitted-nuisance or TMLE estimates.}",
  "\\label{tab:oracle}",
  knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
               digits = c(0, 0, 0, 1, 1, 1),
               col.names = c("Standard error", "Mean", "RMS", "Coverage",
                             "Lower 95\\%", "Upper 95\\%"), escape = FALSE),
  "\\end{table}"
), "paper/tables/semisynth_oracle.tex")
number <- function(x, digits = 1) {
  format(round(x, digits), big.mark = ",", trim = TRUE, scientific = FALSE)
}
values <- c(
  OracleReplicates = number(nrow(draws), 0),
  OracleBias = number(empirical$bias),
  OracleBiasMCSE = number(empirical$bias_mcse),
  OracleSD = number(empirical$sd),
  OracleEfficiencySE = number(empirical$bound_se),
  OracleCoverage = number(100 * empirical$coverage),
  OracleConditionalCoverage = number(100 * conditional$coverage)
)
writeLines(sprintf("\\newcommand{\\%s}{%s}", names(values), values),
           "paper/tables/semisynth_oracle_macros.tex")
