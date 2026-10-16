# Rscript R/recode_report.R
# Summarises data/recode_predictions.csv (from R/recode_collect.R): how much
# test-set predictions move across codings, for HAL and symmetric HAL.
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
dir.create("paper/tables", showWarnings = FALSE)

pred <- read_csv("data/recode_predictions.csv", show_col_types = FALSE) |>
  mutate(method = recode(method, hal = "HAL", shal = "symmetric HAL")) |>
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
