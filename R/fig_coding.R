# Rscript R/fig_coding.R
# One data set, four sign codings of (x1, x2). HAL and symmetric HAL are fit
# under each coding with the same CV folds and every fit is drawn on the
# original axes as its error, fitted minus true value: a coding-invariant
# estimator gives four identical panels.
suppressMessages({
  library(tidyverse)
  library(Matrix)
})
source(file.path("R", "basis.R"))

set.seed(7)
n <- 500
truth <- function(x) 2 * (x[, 1] <= -0.2 & x[, 2] <= -0.2)
x <- matrix(runif(2 * n, -1, 1), n, 2, dimnames = list(NULL, c("x1", "x2")))
y <- truth(x) + rnorm(n, sd = 1.2)
foldid <- sample(rep_len(1:10, n))
grid <- as.matrix(expand_grid(x1 = seq(-1, 1, length.out = 121),
                              x2 = seq(-1, 1, length.out = 121)))

codings <- list("original" = c(1, 1), "x1 negated" = c(-1, 1),
                "x2 negated" = c(1, -1), "both negated" = c(-1, -1))

fit_on_grid <- function(signs, symmetric) {
  xs <- sweep(x, 2, signs, `*`)
  f <- make_dictionary(xs, 19, symmetric)
  cv <- glmnet::cv.glmnet(design_matrix(xs, f), y, standardize = FALSE,
                          foldid = foldid, control = list(thresh = 1e-12))
  as.vector(predict(cv, design_matrix(sweep(grid, 2, signs, `*`), f),
                    s = "lambda.min"))
}

fits <- expand_grid(coding = names(codings), symmetric = c(FALSE, TRUE)) |>
  mutate(fit = map2(coding, symmetric, \(k, s) {
    mutate(as_tibble(grid), fhat = fit_on_grid(codings[[k]], s))
  })) |>
  mutate(method = if_else(symmetric, "symmetric HAL", "HAL"),
         coding = fct_relevel(coding, names(codings))) |>
  unnest(fit) |>
  mutate(f = truth(cbind(x1, x2)))

risk <- fits |>
  summarise(label = sprintf("error %.3f", mean((fhat - f)^2)),
            .by = c(method, coding))
print(risk)

p <- ggplot(fits, aes(x1, x2)) +
  geom_raster(aes(fill = fhat - f)) +
  annotate("rect", xmin = -1, xmax = -0.2, ymin = -1, ymax = -0.2,
           fill = NA, colour = "black", linetype = 2) +
  geom_label(data = risk, aes(label = label), x = 0.97, y = 0.93,
             hjust = 1, size = 2.6, linewidth = 0) +
  facet_grid(method ~ coding) +
  scale_fill_distiller("fitted minus\ntrue value", palette = "RdBu",
                       limits = c(-1, 1), oob = scales::squish) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  coord_equal(expand = FALSE) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme_bw(base_size = 9) +
  theme(panel.spacing = unit(0.6, "lines"))
for (ext in c("pdf", "png")) {
  ggsave(file.path("figures", paste0("coding.", ext)), p, width = 8,
         height = 4.3, dpi = 200)
}
