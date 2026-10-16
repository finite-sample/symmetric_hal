# Rscript R/lalonde_report.R
# Compares the observational estimates in data/lalonde_att.csv (from
# R/lalonde_collect.R) with the experimental benchmark computed from
# data/lalonde_dw.csv, across covariate codings.
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
dir.create("figures", showWarnings = FALSE)
dir.create("paper/tables", showWarnings = FALSE)

method_labels <- c(hal9001 = "HAL, hal9001 defaults",
                   shal9001 = "symmetric HAL, first order",
                   hal = "HAL, zero order", shal = "symmetric HAL, zero order")
experiment <- read_csv("data/lalonde_dw.csv", show_col_types = FALSE) |>
  filter(group %in% c("treated", "experimental")) |>
  summarise(att = mean(re78[treat == 1]) - mean(re78[treat == 0]),
            se = sqrt(var(re78[treat == 1]) / sum(treat == 1) +
                        var(re78[treat == 0]) / sum(treat == 0)))
cat(sprintf("experimental benchmark: %.0f (SE %.0f)\n", experiment$att,
            experiment$se))

# The random sign patterns can repeat (two of 20 did); count each coding once.
att <- read_csv("data/lalonde_att.csv", show_col_types = FALSE) |>
  distinct(method, flipped, .keep_all = TRUE) |>
  mutate(method = factor(method_labels[method], method_labels))
x_range <- c(-8000, 8000)

summary_table <- att |>
  summarise(codings = n(), min = min(att), median = median(att),
            max = max(att), range_se = diff(range(att)) / median(se),
            median_se = median(se),
            median_miss = median(abs(att - experiment$att)), .by = method) |>
  arrange(method)
print(knitr::kable(summary_table, digits = 2))
kbl(summary_table, format = "latex", booktabs = TRUE,
    digits = c(0, 0, 0, 0, 0, 2, 0, 0), format.args = list(big.mark = ","),
    linesep = "", label = "lalonde_att",
    col.names = c("Nuisance fits", "Codings", "Min", "Median", "Max",
                  "Range / SE", "Median SE", "Median distance"),
    caption = paste0("Cross-fitted AIPW estimates of the effect of NSW ",
                     "training on 1978 earnings of the treated (dollars), ",
                     "from the 185 trainees and the PSID-1 comparison group, ",
                     "across random sign codings of the eight covariates. ",
                     "Range / SE is the range divided by the median ",
                     "standard error; median distance is from the ",
                     "experimental estimate, ",
                     format(round(experiment$att), big.mark = ","),
                     " (SE ", round(experiment$se), ").")) |>
  kable_styling(latex_options = c("scale_down", "hold_position")) |>
  save_kable("paper/tables/lalonde_att.tex")

forest <- att |>
  mutate(recode = fct_reorder(factor(recode),
                              att * (method == method_labels[["hal9001"]]),
                              .fun = max)) |>
  ggplot(aes(att, recode)) +
  annotate("rect", xmin = experiment$att - 1.96 * experiment$se,
           xmax = experiment$att + 1.96 * experiment$se, ymin = -Inf,
           ymax = Inf, alpha = 0.12, fill = "steelblue") +
  geom_vline(xintercept = experiment$att, colour = "steelblue") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = att - 1.96 * se, xmax = att + 1.96 * se),
                width = 0, colour = "grey60", orientation = "y") +
  geom_point(size = 0.9) +
  geom_text(aes(x = x_range[1], label = sprintf("%.0f", att)),
            data = \(d) filter(d, att < x_range[1]), hjust = 0, size = 2.3,
            nudge_y = 0.45) +
  coord_cartesian(xlim = x_range) +
  facet_wrap(~method, nrow = 1) +
  labs(x = paste("estimated effect on 1978 earnings of the treated ($),",
                 "95% CI; blue: experimental estimate and its 95% CI;",
                 "off-scale estimates labelled"),
       y = "covariate coding") +
  theme_bw(base_size = 9)
for (ext in c("pdf", "png")) {
  ggsave(file.path("figures", paste0("lalonde_att.", ext)), forest,
         width = 10, height = 4.5, dpi = 150)
}
