# Rscript R/recode_report.R
# Summarises data/recode_predictions.csv and data/recode_ate.csv (from
# R/recode_collect.R): how much predictions and the lalonde treatment effect
# move across codings, for HAL and symmetric HAL.
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
dir.create("figures", showWarnings = FALSE)
dir.create("paper/tables", showWarnings = FALSE)

label <- function(d) {
  mutate(d, method = recode(method, hal = "HAL", shal = "symmetric HAL"))
}
ate <- read_csv("data/recode_ate.csv", show_col_types = FALSE) |> label()
pred <- read_csv("data/recode_predictions.csv", show_col_types = FALSE) |>
  label() |>
  mutate(data = recode(data, boston = "Boston", ozone = "Ozone"))

spread <- pred |>
  group_by(data, method, row) |>
  summarise(codings = n(), sd = sd(pred), range = diff(range(pred)),
            y = first(y), .groups = "drop") |>
  group_by(data) |>
  mutate(sd_y = sd(y)) |>
  group_by(data, method) |>
  summarise(codings = min(codings), mean_sd = mean(sd / sd_y),
            max_range = max(range / sd_y), .groups = "drop")
accuracy <- pred |>
  group_by(data, method, recode) |>
  summarise(rmse = sqrt(mean((pred - y)^2)), .groups = "drop") |>
  group_by(data, method) |>
  summarise(rmse_min = min(rmse), rmse_median = median(rmse),
            rmse_max = max(rmse), .groups = "drop")
prediction_table <- left_join(spread, accuracy, by = c("data", "method"))
print(knitr::kable(prediction_table, digits = 4))

effect_table <- ate |>
  group_by(method) |>
  summarise(codings = n(), ate_min = min(ate), ate_median = median(ate),
            ate_max = max(ate), range_in_se = diff(range(ate)) / median(se),
            sd_in_se = sd(ate) / median(se), median_se = median(se))
print(knitr::kable(effect_table, digits = 2))

kbl(prediction_table, format = "latex", booktabs = TRUE, digits = 3,
    linesep = "", label = "recode_predictions",
    col.names = c("Data", "Method", "Codings", "Mean SD", "Max range",
                  "RMSE min", "RMSE median", "RMSE max"),
    caption = paste("Predictions across random sign codings of the",
                    "covariates, on a fixed test set of 30 percent of the",
                    "rows. Mean SD and max",
                    "range are per test point, in units of the outcome's",
                    "standard deviation.")) |>
  kable_styling(latex_options = "hold_position") |>
  save_kable("paper/tables/recode_predictions.tex")

kbl(effect_table, format = "latex", booktabs = TRUE, digits = 2,
    linesep = "", label = "recode_ate",
    col.names = c("Method", "Codings", "Min", "Median", "Max",
                  "Range / SE", "SD / SE", "Median SE"),
    caption = paste("One-step estimates of the average effect of the NSW",
                    "training program on 1978 earnings (dollars) across",
                    "random sign codings of the eight covariates.")) |>
  kable_styling(latex_options = "hold_position") |>
  save_kable("paper/tables/recode_ate.tex")

forest <- ate |>
  mutate(recode = fct_reorder(factor(recode), ate * (method == "HAL"),
                              .fun = max)) |>
  ggplot(aes(ate, recode)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = ate - 1.96 * se, xmax = ate + 1.96 * se),
                width = 0, colour = "grey60", orientation = "y") +
  geom_point() +
  facet_wrap(~method) +
  labs(x = "estimated average treatment effect on 1978 earnings ($), 95% CI",
       y = "covariate coding") +
  theme_bw(base_size = 9)
for (ext in c("pdf", "png")) {
  ggsave(file.path("figures", paste0("recode_ate.", ext)), forest,
         width = 7, height = 4.5, dpi = 150)
}
