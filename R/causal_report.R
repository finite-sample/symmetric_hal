# Rscript R/causal_report.R
source("R/att.R")
source("R/causal.R")
predictions <- read.csv("data/causal_nuisance.csv")
fit_metadata <- read.csv("data/causal_fits.csv")
groups <- split(predictions, interaction(predictions$method, predictions$coding,
                                         drop = TRUE))
bounds <- range(predictions$y)
results <- do.call(rbind, lapply(groups, function(p) {
  estimates <- compare_att(p, outcome_bounds = bounds)
  data.frame(method = p$method[1], coding = p$coding[1], estimates)
}))
readr::write_csv(results, "data/causal_estimates.csv")
final <- results[results$estimator != "TMLE" |
                   results$step == min(results$step, na.rm = TRUE), ]
stopifnot(all(final$valid_targeting))
final$ci_low <- final$att - 1.96 * final$se
final$ci_high <- final$att + 1.96 * final$se
readr::write_csv(final, "data/causal_summary.csv")
compare_method <- function(method) {
  p <- predictions[predictions$method == method, ]
  original <- p[p$coding == "original", ]
  codings <- setdiff(unique(p$coding), "original")
  do.call(rbind, lapply(codings, function(coding) {
    recoded <- p[p$coding == coding, ]
    comparison <- compare_aipw_codings(original, recoded, bounds)
    data.frame(method = method, coding = coding, comparison)
  }))
}
contrasts <- do.call(rbind, lapply(unique(predictions$method), compare_method))
for (i in seq_len(nrow(contrasts))) {
  d <- contrasts[i, ]
  old <- final[final$method == d$method & final$coding == "original" &
                 final$estimator == d$estimator, ]
  new <- final[final$method == d$method & final$coding == d$coding &
                 final$estimator == d$estimator, ]
  stopifnot(nrow(old) == 1, nrow(new) == 1,
            isTRUE(all.equal(d$att_change, new$att - old$att)),
            isTRUE(all.equal(d$se_change, new$se - old$se)))
}
readr::write_csv(contrasts, "data/causal_coding_contrasts.csv")
display <- final[c("coding", "estimator", "att", "se", "ci_low", "ci_high")]
display <- display[order(match(display$coding,
                               c("original", "continuous_negated", "units")),
                         match(display$estimator,
                               c("AIPW raw", "AIPW", "TMLE"))), ]
labels <- c(original = "Original", continuous_negated = "Signs reversed",
            units = "Different units")
display$coding <- labels[display$coding]
display$estimator[display$estimator == "AIPW"] <- "AIPW, calibrated"
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{ATT estimates under information-preserving covariate recodings.",
  "Nuisance fits use first-order hal9001 defaults and common splits.",
  "Raw AIPW caps propensity scores at 0.99. Calibrated AIPW and TMLE share",
  "clipped outcome predictions and intercept-calibrated initial propensities.",
  "Intervals use the empirical influence-function standard error; their",
  "coverage is not established by this example. The different-units fits",
  "generated convergence warnings and are a software-default diagnostic.}",
  "\\label{tab:causal}",
  knitr::kable(display, format = "latex", booktabs = TRUE, row.names = FALSE,
               digits = c(0, 0, 0, 0, 0, 0),
               col.names = c("Coding", "Estimator", "Estimate", "SE",
                             "Lower 95\\%", "Upper 95\\%"), escape = FALSE),
  "\\end{table}"
), "paper/tables/causal.tex")
contrast_display <- contrasts[c("coding", "estimator", "outcome_component",
                                "propensity_component", "att_change",
                                "absolute_bound", "se_change")]
contrast_display$coding <- labels[contrast_display$coding]
contrast_display$estimator[contrast_display$estimator == "AIPW raw"] <- "Raw"
contrast_display$estimator[contrast_display$estimator == "AIPW"] <- "Calibrated"
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Exact AIPW changes relative to the original coding, in dollars.",
  "The $q$ and $g$ terms in equation~\\eqref{eq:aipw-change} change the",
  "propensity first, then the outcome regression; they sum to $\\Delta$ATT.",
  "Reversing that order reallocates the interaction, so these are not unique",
  "attributions. The bound is $B_q+B_g$. Calibrated rows use the effective",
  "prepared predictions. Different-units nuisance fits report convergence",
  "warnings. These calculations do not decompose TMLE movement.}",
  "\\label{tab:causal-contrast}",
  knitr::kable(contrast_display, format = "latex", booktabs = TRUE,
               row.names = FALSE, digits = c(0, 0, 0, 0, 0, 0, 0),
               col.names = c("Coding", "AIPW", "$q$ term", "$g$ term",
                             "$\\Delta$ATT", "Bound", "$\\Delta$SE"),
               escape = FALSE),
  "\\end{table}"
), "paper/tables/causal_coding_contrasts.tex")
sign_raw <- contrasts[contrasts$coding == "continuous_negated" &
                        contrasts$estimator == "AIPW raw", ]
stopifnot(nrow(sign_raw) == 1)
writeLines(c(
  sprintf("\\newcommand{\\AipwSignMovement}{%.0f}", sign_raw$att_change),
  sprintf("\\newcommand{\\AipwSignSEMovement}{%s}",
          format(round(abs(sign_raw$se_change)), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\AipwSignBound}{%s}",
          format(round(sign_raw$absolute_bound), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\AipwSignPositive}{%s}",
          format(round(sign_raw$positive_row_change), big.mark = ",",
                 trim = TRUE)),
  sprintf("\\newcommand{\\AipwSignNegative}{%s}",
          format(round(sign_raw$negative_row_change), big.mark = ",",
                 trim = TRUE)),
  sprintf("\\newcommand{\\AipwSignMaxG}{%.3f}",
          sign_raw$max_propensity_change),
  sprintf("\\newcommand{\\AipwSignMaxQ}{%s}",
          format(round(sign_raw$max_outcome_change), big.mark = ",",
                 trim = TRUE))
), "paper/tables/causal_coding_macros.tex")
original <- subset(final, coding == "original" & estimator == "TMLE")
reflected <- subset(final, coding == "continuous_negated" & estimator == "TMLE")
stopifnot(nrow(original) == 1, nrow(reflected) == 1)
writeLines(c(
  sprintf("\\newcommand{\\HalVersion}{%s}",
          unique(fit_metadata$hal9001_version)),
  sprintf("\\newcommand{\\TmleVersion}{%s}",
          as.character(utils::packageVersion("tmle"))),
  sprintf("\\newcommand{\\CausalN}{%s}",
          format(length(unique(predictions$row)), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\CausalTreated}{%s}",
          sum(predictions$a[predictions$coding == "original"])),
  sprintf("\\newcommand{\\TmleOriginal}{%s}",
          format(round(original$att), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\TmleReflected}{%s}",
          format(round(reflected$att), big.mark = ",", trim = TRUE)),
  sprintf("\\newcommand{\\TmleMovement}{%s}",
          format(round(abs(reflected$att - original$att)), big.mark = ",",
                 trim = TRUE))
), "paper/tables/causal_macros.tex")
writeLines(trimws(c(
  capture.output(tools::md5sum(c("R/att.R", "R/causal.R", "R/causal_report.R",
                                 "data/causal_nuisance.csv"))),
  capture.output(sessionInfo())
), which = "right"), "data/causal_targeting_session.txt")
path <- read.csv("data/causal_path_check.csv")
stopifnot(nrow(path) == 1, path$master_jerr == 0,
          path$original_inner_jerr < 0, path$extended_inner_jerr == 0,
          path$clamp_movement == 0,
          path$selected_lambda < path$original_inner_last_lambda,
          path$selected_lambda > path$extended_inner_last_lambda)
scientific <- function(x) {
  parts <- strsplit(sprintf("%.2e", x), "e", fixed = TRUE)[[1]]
  sprintf("%s\\times10^{%d}", parts[1], as.integer(parts[2]))
}
writeLines(c(
  sprintf("\\newcommand{\\PathSelectedLambda}{%s}",
          scientific(path$selected_lambda)),
  sprintf("\\newcommand{\\PathLastLambda}{%s}",
          scientific(path$original_inner_last_lambda)),
  sprintf("\\newcommand{\\PathPredictionMovement}{%.0f}",
          path$extended_prediction_movement),
  sprintf("\\newcommand{\\PathMasterMovement}{%s}",
          scientific(path$master_prediction_movement)),
  sprintf("\\newcommand{\\PathMasterKKT}{%s}",
          scientific(path$master_normalized_kkt))
), "paper/tables/causal_path_macros.tex")
print(final[c("coding", "estimator", "att", "se", "targeting_residual_se")],
      row.names = FALSE)
