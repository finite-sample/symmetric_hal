# Rscript R/sim_report.R <name>
# Reads data/<name>_runs.csv (from R/sim_collect.R) and writes the paper's
# simulation figures to figures/ and its ratio table to paper/tables/. The
# other comparisons quoted in the text are printed.
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
tag <- commandArgs(trailingOnly = TRUE)[1]
dir.create("figures", showWarnings = FALSE)
dir.create("paper/tables", showWarnings = FALSE)

method_labels <- c(cart = "CART", hal = "HAL", hal_relaxed = "relaxed HAL",
                   shal = "symmetric HAL",
                   shal_relaxed = "relaxed symmetric HAL", rf = "random forest")
truth_labels <- c(tree = "tree", additive = "additive steps", xor = "XOR",
                  box_upper = "upper box", box_mixed = "mixed box",
                  box_lower = "lower box", box_interior = "interior box",
                  many_boxes = "eight boxes", friedman1 = "Friedman")
runs <- read_csv(file.path("data", paste0(tag, "_runs.csv")),
                 show_col_types = FALSE)
cat("replicates per cell:",
    paste(range(count(runs, truth, snr, method)$n), collapse = "-"), "\n")

save_plot <- function(p, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path("figures", paste0(tag, "_", name, ".", ext)), p,
           width = width, height = height, dpi = 150)
  }
}
truth_facets <- facet_wrap(~truth, scales = "free_y",
                           labeller = as_labeller(truth_labels))

curves <- runs |>
  mutate(truth = factor(truth, names(truth_labels)),
         method = factor(method_labels[method], method_labels)) |>
  summarise(med = median(rel_risk), lo = quantile(rel_risk, 0.25),
            hi = quantile(rel_risk, 0.75), .by = c(truth, snr, method)) |>
  ggplot(aes(snr, med, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08, colour = NA) +
  geom_line() +
  geom_point(size = 0.8) +
  truth_facets +
  scale_x_log10() +
  labs(x = "signal-to-noise ratio", y = "risk / risk of the constant",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")
save_plot(curves, "risk_vs_snr", 9, 8)

# Methods are fit to the same data in each replicate, so compare them by the
# within-replicate log risk ratio.
paired <- function(a, b) {
  runs |>
    filter(method %in% c(a, b)) |>
    select(truth, snr, rep, method, rel_risk) |>
    pivot_wider(names_from = method, values_from = rel_risk) |>
    mutate(truth = factor(truth, names(truth_labels)),
           log_ratio = log(.data[[a]] / .data[[b]]))
}
median_ratio <- function(d) {
  d |>
    summarise(ratio = exp(median(log_ratio)), .by = c(truth, snr)) |>
    pivot_wider(names_from = snr, values_from = ratio) |>
    arrange(truth) |>
    mutate(truth = truth_labels[as.character(truth)]) |>
    rename(`true function` = truth)
}

orientation <- paired("shal", "hal")
p <- ggplot(orientation, aes(factor(snr), log_ratio)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(width = 0.2, size = 0.4, alpha = 0.4) +
  stat_summary(fun = median, geom = "point", colour = "red", size = 1.5) +
  facet_wrap(~truth, labeller = as_labeller(truth_labels)) +
  labs(x = "signal-to-noise ratio",
       y = paste("log(symmetric HAL risk / HAL risk);",
                 "below 0: symmetric HAL better")) +
  theme_bw(base_size = 9)
save_plot(p, "p1_orientation", 9, 6)

orientation |>
  summarise(symmetric_wins = mean(log_ratio < 0), .by = c(truth, snr)) |>
  pivot_wider(names_from = snr, values_from = symmetric_wins) |>
  knitr::kable(digits = 2, caption = "share of replicates shal wins") |>
  print()
for (cmp in list(c("shal_relaxed", "hal_relaxed"), c("cart", "shal_relaxed"))) {
  paired(cmp[1], cmp[2]) |>
    median_ratio() |>
    knitr::kable(digits = 2, caption = paste(cmp, collapse = " / ")) |>
    print()
}
runs |>
  summarise(seconds = median(seconds), .by = method) |>
  knitr::kable(digits = 1) |>
  print()

ratio <- median_ratio(orientation)
print(knitr::kable(ratio, digits = 2))
kbl(ratio, format = "latex", booktabs = TRUE, digits = 2, linesep = "",
    label = paste0(tag, "_ratio"),
    caption = paste("Risk of symmetric HAL relative to HAL, by true function",
                    "(rows) and signal-to-noise ratio (columns): the median",
                    "over 50 replicates of the within-replicate ratio. Below",
                    "1, symmetric HAL is better.")) |>
  kable_styling(latex_options = "hold_position") |>
  add_header_above(c(" " = 1, "signal-to-noise ratio" = ncol(ratio) - 1)) |>
  save_kable(file.path("paper/tables", paste0(tag, "_ratio.tex")))
