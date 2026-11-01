# Rscript R/semisynth_paper.R
estimates <- read.csv("data/semisynth_att.csv")
methods <- c("oracle", "hal", "canonical", "shal", "hal9001",
             "hal9001_canonical")
if (!setequal(estimates$method, methods) ||
      !setequal(estimates$replicate, 1:20)) {
  stop("paper tables require all six methods and pilot replicates 1 through 20")
}
source(file.path("R", "semisynth_report.R"), local = TRUE)

method_labels <- c("Oracle", "Anchored (0th)", "Canonical (0th)",
                   "Symmetric (0th)", "HAL package (1st)", "Canonical (1st)")
names(method_labels) <- methods
estimators <- c("AIPW raw", "AIPW", "TMLE")
estimator_labels <- c("Raw AIPW", "Prepared AIPW", "TMLE")
names(estimator_labels) <- estimators
codings <- c("original", "continuous_negated")
coding_labels <- c("Original", "Reversed")
names(coding_labels) <- codings
directory <- file.path("paper", "tables")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)

uncertainty <- function(value, error) {
  ifelse(is.na(error), sprintf("%.2f (--)", value),
         sprintf("%.2f (%.2f)", value, error))
}
small_number <- function(x) {
  vapply(x, function(value) {
    if (value == 0) return("0")
    if (abs(value) >= 0.01 && abs(value) < 1e5) {
      return(sprintf("%.2f", value))
    }
    parts <- strsplit(sprintf("%.2e", value), "e", fixed = TRUE)[[1]]
    sprintf("$%s\\times10^{%d}$", parts[1], as.integer(parts[2]))
  }, character(1))
}
write_table <- function(display, caption, name, columns,
                        groups = display[[1]], placement = "p") {
  linesep <- ifelse(head(groups, -1) != tail(groups, -1),
                    "\\addlinespace", "")
  writeLines(c(
    paste0("\\begin{table}[", placement, "]"), "\\centering", "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:semisynth-", name, "}"),
    knitr::kable(display, format = "latex", booktabs = TRUE,
                 row.names = FALSE, col.names = columns, escape = FALSE,
                 linesep = linesep,
                 align = c("l", "l", rep("r", ncol(display) - 2))),
    "\\end{table}"
  ), file.path(directory, paste0("semisynth_", name, ".tex")))
}

for (estimator in estimators) {
  data <- summary[summary$estimator == estimator, ]
  row_order <- order(match(data$method, methods), match(data$coding, codings))
  data <- data[row_order, ]
  display <- data.frame(
    learner = unname(method_labels[data$method]),
    coding = unname(coding_labels[data$coding]),
    bias = uncertainty(data$bias / 1000, data$bias_mcse / 1000),
    rmse = uncertainty(data$rmse / 1000, data$rmse_mcse / 1000),
    coverage = sprintf("%.0f [%.0f, %.0f]", 100 * data$coverage,
                       100 * data$coverage_lower, 100 * data$coverage_upper),
    warnings = data$nuisance_warning_replicates
  )
  caption <- paste0(
    estimator_labels[estimator], " accuracy in the known-effect pilot. ",
    "Each row uses the same 20 draws of ", unique(estimates$n),
    " observations from the fixed empirical-covariate population; the ATT is ",
    unique(estimates$truth), " dollars. Bias and RMSE are in thousands of ",
    "dollars, with Monte Carlo standard errors in parentheses (delta method ",
    "for RMSE). Coverage is the percentage of nominal 95\\% intervals ",
    "containing the ATT, with an exact 95\\% binomial interval in brackets. ",
    "Warnings count draws with any nuisance-fit warning; these draws are ",
    "retained. Oracle denotes true initial nuisances; preparation changes ",
    "them for prepared AIPW and TMLE. Reversed negates continuous covariates ",
    "using the same draws and folds. 0th and 1st denote basis order. ",
    "Twenty draws do not establish calibrated coverage. A missing RMSE ",
    "MCSE (--) indicates zero empirical RMSE."
  )
  name <- c("aipw_raw", "aipw", "tmle")[match(estimator, estimators)]
  write_table(display, caption, name,
              c("Learner", "Coding", "Bias (MCSE)", "RMSE (MCSE)",
                "Coverage [CI]", "Warnings"),
              placement = if (name == "aipw_raw") "htbp" else "p")
}

data <- movement[order(match(movement$method, methods),
                       match(movement$estimator, estimators)), ]
display <- data.frame(
  learner = unname(method_labels[data$method]),
  estimator = unname(estimator_labels[data$estimator]),
  mean = paste0(small_number(data$mean_movement), " (",
                small_number(data$movement_mcse), ")"),
  maximum = small_number(data$maximum_movement),
  se_maximum = small_number(data$maximum_se_movement)
)
write_table(display, paste0(
  "Coding movement in the known-effect pilot. Each learner-estimator pair ",
  "uses all 20 matched draws under original and reversed continuous ",
  "covariates, with common folds. Entries are absolute differences in ",
  "dollars. The mean ATT movement has its Monte Carlo standard error in ",
  "parentheses; maxima describe these draws and have no uncertainty interval. ",
  "SE movement concerns empirical influence-function standard errors. ",
  "Small nonzero numerical differences are retained. Agreement under ",
  "recoding does not establish accuracy or correct interval coverage."
), "movement", c("Learner", "Estimator", "Mean ATT (MCSE)",
                 "Max ATT", "Max SE"))

pair_names <- c("1st: canonical / package", "0th: canonical / anchored",
                "0th: symmetric / anchored", "0th: symmetric / canonical")
pair_keys <- paste(comparators$method, comparators$reference)
data <- comparisons
index <- match(paste(data$method, data$reference), pair_keys)
data$comparison <- pair_names[index]
data$ratio_interval <- ifelse(
  is.na(data$rmse_ratio), "--",
  ifelse(is.na(data$rmse_ratio_lower), sprintf("%.2f [--]", data$rmse_ratio),
         sprintf("%.2f [%.2f, %.2f]", data$rmse_ratio,
                 data$rmse_ratio_lower, data$rmse_ratio_upper))
)
original <- data[data$coding == "original", ]
reversed <- data[data$coding == "continuous_negated", ]
paired <- merge(original, reversed, by = c("comparison", "estimator"))
paired <- paired[order(match(paired$estimator, estimators),
                       match(paired$comparison, pair_names)), ]
display <- data.frame(
  comparison = paired$comparison,
  estimator = unname(estimator_labels[paired$estimator]),
  original = paired$ratio_interval.x, reversed = paired$ratio_interval.y
)
write_table(display, paste0(
  "Paired accuracy comparisons in the known-effect pilot. Each ratio is ",
  "the RMSE of the numerator learner divided by that of the denominator, ",
  "using all 20 common draws. Values below one favor the numerator. ",
  "Brackets give pointwise approximate 95\\% Monte Carlo intervals from ",
  "paired squared errors, using a delta-method log-ratio standard error ",
  "and a $t_{19}$ critical value. These exploratory intervals are not ",
  "multiplicity-adjusted and do not describe uncertainty about the ATT. ",
  "An unavailable ratio (--) indicates zero reference RMSE; an unavailable ",
  "interval [--] indicates zero numerator RMSE."
), "comparisons", c("RMSE numerator / denominator", "Estimator",
                    "Original [CI]", "Reversed [CI]"),
groups = display$estimator)

default <- subset(movement, method == "hal9001" & estimator == "TMLE")
symmetric <- subset(movement, method == "shal" & estimator == "TMLE")
oracle <- subset(summary, method == "oracle" & coding == "original" &
                   estimator == "AIPW raw")
ratio_interval <- function(method, reference, coding) {
  row <- comparisons[comparisons$method == method &
                       comparisons$reference == reference &
                       comparisons$coding == coding &
                       comparisons$estimator == "TMLE", ]
  stopifnot(nrow(row) == 1)
  sprintf("%.2f [%.2f, %.2f]", row$rmse_ratio,
          row$rmse_ratio_lower, row$rmse_ratio_upper)
}
macros <- c(
  CausalPilotDraws = length(unique(estimates$replicate)),
  CausalPilotDefaultMovement = sprintf("%.0f", default$mean_movement),
  CausalPilotDefaultMovementMCSE = sprintf("%.0f", default$movement_mcse),
  CausalPilotSymmetricMaximum = small_number(symmetric$maximum_movement),
  CausalPilotFirstOriginal = ratio_interval("hal9001_canonical", "hal9001",
                                            "original"),
  CausalPilotFirstReversed = ratio_interval("hal9001_canonical", "hal9001",
                                            "continuous_negated"),
  CausalPilotSymOriginal = ratio_interval("shal", "hal", "original"),
  CausalPilotSymReversed = ratio_interval("shal", "hal", "continuous_negated"),
  CausalPilotOracleBias = sprintf("%.0f", oracle$bias),
  CausalPilotOracleMCSE = sprintf("%.0f", oracle$bias_mcse)
)
writeLines(sprintf("\\newcommand{\\%s}{%s}", names(macros), macros),
           file.path(directory, "semisynth_macros.tex"))
