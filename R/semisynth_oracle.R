# Rscript R/semisynth_oracle.R
source("R/att.R")
source("R/semisynth.R")
dgp <- semisynth_dgp(read.csv("data/lalonde_dw.csv"))
results <- lapply(20001:30000, function(replicate) {
  sim <- semisynth_draw(dgp, replicate)
  g <- sim$truth$true_g
  fit <- aipw_att(sim$a, sim$y, sim$truth$true_q0, g)
  conditional_se <- dgp$sigma * sqrt(sum(
    sim$a + (1 - sim$a) * (g / (1 - g))^2
  )) / sum(sim$a)
  data.frame(replicate = replicate, att = fit$att, se = fit$se,
             conditional_se = conditional_se, truth = sim$effect,
             bound_se = sim$bound_se, treated = sum(sim$a))
})
estimates <- do.call(rbind, results)
readr::write_csv(estimates, "data/semisynth_oracle.csv")
error <- estimates$att - estimates$truth
summary <- do.call(rbind, lapply(c("se", "conditional_se"), function(column) {
  covered <- abs(error) <= stats::qnorm(0.975) * estimates[[column]]
  interval <- stats::binom.test(sum(covered), nrow(estimates))$conf.int
  data.frame(interval = column, replicates = nrow(estimates),
             bias = mean(error), bias_mcse = sd(error) / sqrt(nrow(estimates)),
             sd = sd(estimates$att), rmse = sqrt(mean(error^2)),
             mean_se = mean(estimates[[column]]),
             rms_se = sqrt(mean(estimates[[column]]^2)),
             coverage = mean(covered),
             coverage_lower = interval[1], coverage_upper = interval[2],
             bound_se = sqrt(dgp$efficiency_variance / nrow(dgp$population)))
}))
readr::write_csv(summary, "data/semisynth_oracle_summary.csv")
writeLines(trimws(c(capture.output(tools::md5sum(c(
  "R/att.R", "R/semisynth.R", "R/semisynth_oracle.R", "data/lalonde_dw.csv"
))), capture.output(sessionInfo())), which = "right"),
"data/semisynth_oracle_session.txt")
print(summary, row.names = FALSE)
