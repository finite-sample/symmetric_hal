# Rscript R/audit_report.R
# How far each estimator's test-set predictions move under each recoding of
# Boston housing (data/audit_predictions.csv, from R/audit_collect.R), in
# standard deviations of the outcome.
suppressMessages({
  library(tidyverse)
  library(kableExtra)
})
dir.create("paper/tables", showWarnings = FALSE)

families <- tribble(
  ~method,              ~label,                      ~family,
  "ols",                "OLS",                       "Unpenalized",
  "ols_interactions",   "OLS, pairwise interactions", "Unpenalized",
  "gam",                "GAM (mgcv)",                "Unpenalized",
  "ridge",              "ridge",                     "Penalized (glmnet)",
  "lasso",              "lasso",                     "Penalized (glmnet)",
  "elastic_net",        "elastic net",               "Penalized (glmnet)",
  "lasso_interactions", "lasso, pairwise interactions", "Penalized (glmnet)",
  "cart",               "CART (rpart)",              "Trees",
  "random_forest",      "random forest (ranger)",    "Trees",
  "boosting",           "boosting (xgboost)",        "Trees",
  "bart",               "BART (dbarts)",             "Trees",
  "knn",                "10 nearest neighbours",     "Other",
  "svm",                "SVM (e1071)",               "Other",
  "mars",               "MARS (earth)",              "Other",
  "hal",                "HAL, zero order",           "HAL",
  "shal",               "symmetric HAL, zero order", "HAL",
  "hal9001",            "HAL, hal9001 defaults",     "HAL",
  "shal9001",           "symmetric HAL, first order", "HAL",
  "shal9001_scaled",    "symmetric HAL, first order, scaled", "HAL"
)
codings <- c(negated = "Negated", negated_from_max = "Max minus x",
             shifted = "Shifted", rescaled = "Units",
             reference_level = "Reference level", reseed = "New seed")

pred <- read_csv("data/audit_predictions.csv", show_col_types = FALSE)
sd_y <- sd(filter(pred, method == "ols", coding == "original")$y)
moves <- pred |>
  filter(coding != "original") |>
  left_join(filter(pred, coding == "original") |>
              select(method, row, original = pred),
            by = c("method", "row")) |>
  summarise(median = median(abs(pred - original)) / sd_y,
            max = max(abs(pred - original)) / sd_y,
            .by = c(method, coding))
print(knitr::kable(moves |> arrange(method, coding), digits = 3))

# Below 0.005 SD at every test point the change is solver tolerance.
cell <- function(median, max) {
  if_else(max < 0.005, "0", sprintf("%.2f (%.2f)", median, max))
}
rmse <- pred |>
  filter(coding == "original") |>
  summarise(rmse = sqrt(mean((pred - y)^2)), .by = method)
table <- moves |>
  mutate(cell = cell(median, max), coding = codings[coding]) |>
  select(method, coding, cell) |>
  pivot_wider(names_from = coding, values_from = cell) |>
  left_join(rmse, by = "method")
table <- families |>
  left_join(table, by = "method") |>
  select(label, all_of(unname(codings)), rmse, family)

kbl(select(table, -family), format = "latex", booktabs = TRUE, digits = 2,
    linesep = "", label = "audit", align = c("l", rep("r", 7)),
    col.names = c("Estimator", unname(codings), "Test RMSE"),
    caption = paste("Boston housing: median (and largest) change in a",
                    "test-set prediction when the covariates are recoded,",
                    "in standard deviations of the outcome; 0 means under",
                    "0.005 at every test point. Negated, max minus x,",
                    "shifted (by ten standard deviations) and units",
                    "(times 100, 0.01 or 12) recode all eleven continuous",
                    "covariates; reference level changes the omitted",
                    "category of the nine-level highway-access index; new",
                    "seed refits the original coding with another random",
                    "seed. Test RMSE is for the original coding.")) |>
  pack_rows(index = table(fct_inorder(table$family))) |>
  kable_styling(latex_options = c("scale_down", "hold_position")) |>
  save_kable("paper/tables/audit.tex")
