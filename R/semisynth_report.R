# Rscript R/semisynth_report.R
suppressPackageStartupMessages(library(dplyr))
estimates <- read.csv("data/semisynth_att.csv")
nuisance_metrics <- c("q0_risk", "q0_treated_risk", "q1_treated_risk", "g_risk",
                      "raw_drift", "prepared_drift")
clipping <- c("prepared_clipped_q0", "prepared_clipped_q1",
              "prepared_clipped_g")
job_metrics <- c(nuisance_metrics, clipping, "nuisance_warnings",
                 "nuisance_seconds", "targeting_seconds", "seconds", "n",
                 "truth", "bound_se", "q_targeting_floor")
stopifnot(all(c(job_metrics, "att", "se", "valid_targeting") %in%
                names(estimates)),
          all(is.finite(as.matrix(estimates[c(job_metrics, "att", "se")]))),
          all(estimates$se >= 0), all(estimates$bound_se > 0),
          !anyDuplicated(estimates[c("replicate", "method", "coding",
                                     "estimator")]),
          all(estimates$valid_targeting),
          length(unique(estimates$truth)) == 1,
          length(unique(estimates$bound_se)) == 1,
          length(unique(estimates$n)) == 1,
          length(unique(estimates$q_targeting_floor)) == 1,
          setequal(estimates$coding, c("original", "continuous_negated")),
          setequal(estimates$estimator, c("AIPW raw", "AIPW", "TMLE")))
replicates <- sort(unique(estimates$replicate))
stopifnot(length(replicates) >= 2,
          nrow(estimates) == length(replicates) *
            length(unique(estimates$method)) * 2 * 3)
job_columns <- c("replicate", "method", "coding", job_metrics)
jobs <- unique(estimates[job_columns])
stopifnot(!anyDuplicated(jobs[c("replicate", "method", "coding")]))
groups <- split(estimates, interaction(estimates$method, estimates$coding,
                                       estimates$estimator, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(data) {
  error <- data$att - data$truth
  covered <- abs(error) <= 1.96 * data$se
  count <- nrow(data)
  interval <- stats::binom.test(sum(covered), count)$conf.int
  rmse <- sqrt(mean(error^2))
  data.frame(data[1, c("method", "coding", "estimator")],
             replicates = count, bias = mean(error),
             bias_mcse = sd(error) / sqrt(count), sd = sd(data$att),
             rmse = rmse, rmse_mcse = if (rmse > 0) {
               sd(error^2) / (2 * rmse * sqrt(count))
             } else {
               NA_real_
             },
             mean_se = mean(data$se), rms_se = sqrt(mean(data$se^2)),
             bound_se = unique(data$bound_se),
             coverage = mean(covered), coverage_lower = interval[1],
             coverage_upper = interval[2],
             nuisance_warning_replicates = sum(data$nuisance_warnings > 0),
             row.names = NULL)
}))
readr::write_csv(summary, "data/semisynth_summary.csv")
original <- subset(estimates, coding == "original")
reflected <- subset(estimates, coding == "continuous_negated")
paired <- merge(original, reflected, by = c("replicate", "method", "estimator"))
stopifnot(nrow(paired) == nrow(original), nrow(paired) == nrow(reflected))
paired$movement <- abs(paired$att.x - paired$att.y)
paired$se_movement <- abs(paired$se.x - paired$se.y)
movement <- paired |>
  group_by(method, estimator) |>
  summarise(replicates = n(), mean_movement = mean(movement),
            movement_mcse = sd(movement) / sqrt(n()),
            median_movement = median(movement),
            maximum_movement = max(movement),
            mean_se_movement = mean(se_movement),
            maximum_se_movement = max(se_movement),
            .groups = "drop")
readr::write_csv(movement, "data/semisynth_movement.csv")
nuisance <- jobs |>
  group_by(method, coding) |>
  summarise(across(all_of(nuisance_metrics),
                   list(mean = mean,
                        mcse = function(z) sd(z) / sqrt(length(z)))),
            .groups = "drop")
readr::write_csv(nuisance, "data/semisynth_nuisance_summary.csv")
diagnostics <- jobs |>
  group_by(method, coding) |>
  summarise(replicates = n(),
            across(all_of(clipping),
                   list(mean = mean, replicates = function(z) sum(z > 0))),
            nuisance_warning_replicates = sum(nuisance_warnings > 0),
            mean_nuisance_seconds = mean(nuisance_seconds),
            mean_targeting_seconds = mean(targeting_seconds),
            mean_seconds = mean(seconds), .groups = "drop")
readr::write_csv(diagnostics, "data/semisynth_diagnostics.csv")

comparators <- data.frame(
  method = c("hal9001_canonical", "canonical", "shal", "shal"),
  reference = c("hal9001", "hal", "hal", "canonical")
)
comparisons <- list()
for (i in seq_len(nrow(comparators))) {
  method <- comparators$method[i]
  reference <- comparators$reference[i]
  if (!all(c(method, reference) %in% estimates$method)) next
  matched <- merge(estimates[estimates$method == method, ],
                   estimates[estimates$method == reference, ],
                   by = c("replicate", "coding", "estimator"))
  cells <- split(matched, interaction(matched$coding, matched$estimator,
                                      drop = TRUE))
  compare_cell <- function(p) {
    error <- p$att.x - p$truth.x
    reference_error <- p$att.y - p$truth.y
    stopifnot(identical(sort(p$replicate), replicates),
              identical(p$truth.x, p$truth.y))
    count <- nrow(p)
    mse <- mean(error^2)
    reference_mse <- mean(reference_error^2)
    difference <- error^2 - reference_error^2
    difference_mcse <- sd(difference) / sqrt(count)
    critical <- stats::qt(0.975, count - 1)
    log_ratio <- ratio_mcse <- NA_real_
    ratio <- if (reference_mse > 0) sqrt(mse / reference_mse) else NA_real_
    if (mse > 0 && reference_mse > 0) {
      log_ratio <- 0.5 * log(mse / reference_mse)
      ratio_mcse <- 0.5 * sd(error^2 / mse -
                               reference_error^2 / reference_mse) / sqrt(count)
    }
    data.frame(
      method = method, reference = reference,
      coding = p$coding[1], estimator = p$estimator[1], replicates = count,
      mse = mse, reference_mse = reference_mse,
      mse_difference = mean(difference), mse_difference_mcse = difference_mcse,
      mse_difference_lower = mean(difference) - critical * difference_mcse,
      mse_difference_upper = mean(difference) + critical * difference_mcse,
      rmse_ratio = ratio, log_rmse_ratio_mcse = ratio_mcse,
      rmse_ratio_lower = exp(log_ratio - critical * ratio_mcse),
      rmse_ratio_upper = exp(log_ratio + critical * ratio_mcse)
    )
  }
  rows <- lapply(cells, compare_cell)
  comparisons[[length(comparisons) + 1]] <- do.call(rbind, rows)
}
empty_comparisons <- data.frame(
  method = character(), reference = character(), coding = character(),
  estimator = character(), replicates = integer(), mse = numeric(),
  reference_mse = numeric(), mse_difference = numeric(),
  mse_difference_mcse = numeric(), mse_difference_lower = numeric(),
  mse_difference_upper = numeric(), rmse_ratio = numeric(),
  log_rmse_ratio_mcse = numeric(), rmse_ratio_lower = numeric(),
  rmse_ratio_upper = numeric()
)
comparisons <- bind_rows(empty_comparisons, comparisons)
readr::write_csv(comparisons, "data/semisynth_comparisons.csv")
print(summary, row.names = FALSE)
print(movement)
