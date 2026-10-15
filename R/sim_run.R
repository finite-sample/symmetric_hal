# Rscript R/sim_run.R <design> <reps> <out_dir> <cores>
# Each replicate is saved as soon as it finishes, so an interrupted run
# resumes where it stopped; a replicate that errors is saved with its error.
suppressMessages(library(Matrix))
for (f in c("basis", "fit", "sim_setup")) {
  source(file.path("R", paste0(f, ".R")))
}

args <- commandArgs(trailingOnly = TRUE)
design <- args[1]
reps <- seq_len(as.integer(args[2]))
out_dir <- args[3]
cores <- as.integer(args[4])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

snrs <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4)
grid <- switch(design,
  pilot = expand.grid(truth = names(truths), snr = 1, n = 500, rho = 0,
                      stringsAsFactors = FALSE),
  tier1 = expand.grid(truth = names(truths), snr = snrs, n = 500, rho = 0,
                      stringsAsFactors = FALSE),
  tier2 = rbind(
    expand.grid(truth = c("tree", "box_lower", "many_boxes", "friedman1"),
                snr = snrs, n = c(200, 2000), rho = 0,
                stringsAsFactors = FALSE),
    expand.grid(truth = c("tree", "box_lower", "many_boxes", "friedman1"),
                snr = snrs, n = 500, rho = 0.5, stringsAsFactors = FALSE)
  )
)
# Replicate-major order, so a run stopped early leaves every cell with the
# same number of finished replicates.
grid <- merge(grid, data.frame(rep = reps))
grid <- grid[order(grid$rep), ]
d <- 10

test_sets <- lapply(split(grid, grid$rho), function(g) {
  set.seed(99)
  draw_x(2000, d, g$rho[1])
})

path_of <- function(g) {
  file.path(out_dir, sprintf("%s_snr%s_n%d_rho%s_rep%03d.rds", g$truth,
                             g$snr, g$n, g$rho, g$rep))
}
todo <- which(!file.exists(vapply(seq_len(nrow(grid)), function(i) {
  path_of(grid[i, ])
}, character(1))))
cat(length(todo), "of", nrow(grid), "replicates to run\n")

invisible(parallel::mclapply(todo, function(i) {
  g <- grid[i, ]
  result <- tryCatch(
    run_replicate(g$truth, g$snr, g$rep, g$n, d, g$rho,
                  test_sets[[as.character(g$rho)]], methods_all),
    error = function(e) list(error = conditionMessage(e), cell = g)
  )
  saveRDS(result, path_of(g))
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
