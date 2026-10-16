# Rscript R/semisynth_report.R
# Bias, variability against the efficiency bound, coverage and movement
# across codings for the semi-synthetic LaLonde estimates in
# data/semisynth_att.csv (from R/semisynth_collect.R).
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
dir.create("figures", showWarnings = FALSE)
dir.create("paper/tables", showWarnings = FALSE)

method_labels <- c(hal9001 = "HAL, hal9001 defaults",
                   shal9001 = "symmetric HAL, first order")
coding_labels <- c(original = "original",
                   continuous_negated = "continuous covariates negated")
est <- read_csv("data/semisynth_att.csv", show_col_types = FALSE)
if ("error" %in% names(est)) {
  cat(sum(!is.na(est$error)), "estimates errored\n")
  est <- filter(est, is.na(error))
}
est <- est |>
  mutate(method = factor(method_labels[method], method_labels),
         coding = factor(coding_labels[coding], coding_labels))

performance <- est |>
  summarise(reps = n(), bias = mean(att - truth), sd = sd(att),
            rmse = sqrt(mean((att - truth)^2)), sd_over_bound = sd /
              mean(bound_se), mean_se = mean(se),
            coverage = mean(abs(att - truth) <= 1.96 * se),
            .by = c(method, coding)) |>
  arrange(method, coding)
movement <- est |>
  select(rep, method, coding, att) |>
  pivot_wider(names_from = coding, values_from = att) |>
  mutate(move = abs(.data[[coding_labels[[1]]]] -
                      .data[[coding_labels[[2]]]])) |>
  summarise(median_move = median(move), .by = method)
table <- left_join(performance, movement, by = "method")
print(knitr::kable(table, digits = 2))

kbl(table, format = "latex", booktabs = TRUE,
    digits = c(0, 0, 0, 0, 0, 0, 2, 0, 2, 0),
    format.args = list(big.mark = ","), linesep = "",
    label = "semisynth",
    col.names = c("Nuisance fits", "Coding", "Reps", "Bias", "SD", "RMSE",
                  "SD / bound", "Mean SE", "Coverage", "Median move"),
    caption = paste("Semi-synthetic LaLonde: estimates of an effect on the",
                    "treated of 1,800 dollars. SD / bound divides the",
                    "standard deviation across replicates by the",
                    "efficiency bound; coverage is that of the nominal 95",
                    "percent interval; median move is the median absolute",
                    "change in the estimate between the two codings of the",
                    "same replicate.")) |>
  kable_styling(latex_options = c("scale_down", "hold_position")) |>
  save_kable("paper/tables/semisynth.tex")

p <- ggplot(est, aes(coding, att)) +
  geom_hline(yintercept = unique(est$truth), colour = "steelblue") +
  geom_line(aes(group = rep), colour = "grey80", linewidth = 0.3) +
  geom_point(size = 0.6, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", colour = "red", size = 2) +
  facet_wrap(~method) +
  labs(x = NULL, y = paste("estimated effect on the treated ($); blue:",
                           "truth, red: mean, lines join a replicate")) +
  theme_bw(base_size = 9)
for (ext in c("pdf", "png")) {
  ggsave(file.path("figures", paste0("semisynth.", ext)), p, width = 8,
         height = 4.5, dpi = 150)
}
