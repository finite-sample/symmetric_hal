root <- file.path("..", "..")
for (f in c("basis", "fit")) {
  source(file.path(root, "R", paste0(f, ".R")))
}

toy_x <- function(n = 60, d = 3, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * d, -1, 1), n, d)
  colnames(x) <- paste0("x", seq_len(d))
  x
}
