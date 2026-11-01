# Rscript R/pilot_report.R [data_dir] [figure_dir] [table_dir]
suppressPackageStartupMessages(library(tidyverse))
source("R/pilot.R")
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
figure_dir <- if (length(args) >= 2) args[2] else "figures"
table_dir <- if (length(args) >= 3) args[3] else "paper/tables"
comparisons <- read_csv(file.path(data_dir, "pilot_comparisons.csv"),
                        show_col_types = FALSE)
summary <- read_csv(file.path(data_dir, "pilot_summary.csv"),
                    show_col_types = FALSE)
records <- read_csv(file.path(data_dir, "pilot_runs.csv"),
                    show_col_types = FALSE)
diagnostics <- read_csv(file.path(data_dir, "pilot_diagnostics.csv"),
                        show_col_types = FALSE)
stopifnot(nrow(diagnostics) == 2 * length(unique(records$job)),
          !anyDuplicated(diagnostics[c("job", "method")]),
          all(is.finite(diagnostics$max_reflection_movement)),
          all(diagnostics$max_reflection_movement >= 0),
          setequal(diagnostics$method, c("symmetric", "canonical")))
for (method in c("symmetric", "canonical")) {
  stopifnot(setequal(diagnostics$job[diagnostics$method == method],
                     records$job))
}
keys <- c("truth", "d", "n", "snr", "distribution", "knot_rule")
expected <- distinct(pilot_grid(1), across(all_of(keys)))
methods <- c("anchored", "selected", "averaged", "canonical", "cart", "rf")
groups <- merge(expected, data.frame(comparator = methods), by = NULL)
stopifnot(nrow(anti_join(groups, comparisons, by = c(keys, "comparator"))) == 0,
          nrow(anti_join(comparisons, groups, by = c(keys, "comparator"))) == 0,
          !anyDuplicated(comparisons[c(keys, "comparator")]),
          all(comparisons$replicates == 20),
          all(is.finite(comparisons$log_ratio_mcse)),
          all(comparisons$ratio_lower > 0),
          all(comparisons$risk_ratio > 0))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
labels <- c(anchored = "Anchored", selected = "Selected anchor",
            averaged = "Anchor average", canonical = "Canonical orientation")
colors <- c(anchored = "#0072B2", selected = "#D55E00",
            averaged = "#009E73", canonical = "#CC79A7")
truth_labels <- c(upper = "Upper orthant", lower = "Lower orthant",
                  mixed = "Mixed orthant", interior = "Interior box",
                  conflicting = "Conflicting orientations",
                  additive = "Additive steps", smooth = "Smooth signal")
primary <- comparisons |>
  filter(comparator %in% names(labels), knot_rule == "quantile") |>
  mutate(comparator = factor(comparator, levels = names(labels)),
         distribution = factor(
           distribution, levels = c("independent", "correlated", "discrete")
         ),
         d = paste0("Dimension ", d), snr = factor(snr))
position <- position_dodge(width = 0.65)
for (truth in names(truth_labels)) {
  data <- filter(primary, .data$truth == .env$truth)
  plot <- ggplot(data, aes(x = factor(n), y = risk_ratio, color = comparator,
                           shape = snr, linetype = snr,
                           group = interaction(comparator, snr))) +
    geom_hline(yintercept = 1, color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ratio_lower, ymax = ratio_upper),
                  width = 0.08, position = position, linewidth = 0.35) +
    geom_line(position = position, linewidth = 0.35) +
    geom_point(position = position, size = 1.8) +
    facet_grid(d ~ distribution) +
    scale_y_log10() +
    scale_color_manual(values = colors, labels = labels, drop = FALSE) +
    labs(x = "Training observations", y = "Mean risk: symmetric / comparator",
         title = truth_labels[[truth]], color = "Comparator", shape = "SNR",
         linetype = "SNR",
         caption = paste(
           "Below 1 favors symmetric HAL. Bars: approximate 95% Monte Carlo",
           "intervals from paired draws.\nQuantile knots;",
           "20 replicates per cell.",
           "Uniform/discrete risk is integrated exactly; correlated risk uses",
           "10,000 common test points."
         )) +
    guides(color = guide_legend(nrow = 2)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 11),
          axis.text = element_text(size = 11),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, size = 10))
  for (extension in c("pdf", "png")) {
    ggsave(file.path(figure_dir, paste0("pilot_", truth, ".", extension)),
           plot, width = 8.5, height = 6.8, dpi = 180)
  }
}
table <- summary |>
  filter(knot_rule == "fixed", method %in% c(names(labels), "symmetric")) |>
  arrange(n, match(method, c(names(labels), "symmetric"))) |>
  select(n, method, risk, risk_mcse, mean_seconds)
table$method <- c(labels, symmetric = "Symmetric")[table$method]
writeLines(c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Prediction risk in the fixed-knot mechanism comparison.",
  "The signal combines two boxes with conflicting preferred orientations;",
  "dimension is two, covariates are independent uniform, and SNR is two.",
  "Each row averages 20 training samples. Risk is integrated squared error;",
  "MCSE measures uncertainty over simulation draws. Fitting time includes",
  "tuning and refitting, excluding the study's additional reflection checks.",
  "Times are elapsed seconds under concurrent machine load, not isolated",
  "serial benchmarks.",
  "All methods use the same knots, penalty candidates, and folds.}",
  "\\label{tab:pilot-fixed}",
  knitr::kable(table, format = "latex", booktabs = TRUE, row.names = FALSE,
               digits = c(0, 0, 5, 5, 2),
               col.names = c("$n$", "Method", "Mean risk", "MCSE", "Seconds"),
               escape = FALSE), "\\end{table}"
), file.path(table_dir, "pilot_fixed.tex"))
fixed <- filter(comparisons, knot_rule == "fixed")
ratio_macros <- unlist(lapply(c(Small = 100, Large = 300), function(n) {
  size <- if (n == 100) "Small" else "Large"
  unlist(lapply(c("selected", "averaged", "canonical"), function(method) {
    name <- c(selected = "Selected", averaged = "Averaged",
              canonical = "Canonical")[[method]]
    row <- filter(fixed, .data$n == .env$n, comparator == method)
    stopifnot(nrow(row) == 1)
    sprintf("\\newcommand{\\Pilot%s%s%s}{%.2f}", name, size,
            c("Ratio", "Lower", "Upper"),
            c(row$risk_ratio, row$ratio_lower, row$ratio_upper))
  }))
}))
scientific <- function(x) {
  parts <- strsplit(sprintf("%.1e", x), "e", fixed = TRUE)[[1]]
  sprintf("%s\\times10^{%d}", parts[1], as.integer(parts[2]))
}
writeLines(c(
  sprintf("\\newcommand{\\PilotCases}{%s}", length(unique(records$job))),
  sprintf("\\newcommand{\\PilotScenarios}{%s}", nrow(expected)),
  sprintf("\\newcommand{\\PilotReplicates}{%s}",
          unique(comparisons$replicates)),
  sprintf("\\newcommand{\\PilotMaxMovement}{%s}", scientific(max(
    diagnostics$max_reflection_movement[diagnostics$method == "symmetric"]
  ))),
  sprintf("\\newcommand{\\PilotCanonicalMovement}{%s}", scientific(max(
    diagnostics$max_reflection_movement[diagnostics$method == "canonical"]
  ))), ratio_macros
), file.path(table_dir, "pilot_macros.tex"))
