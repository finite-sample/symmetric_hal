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

column_keys <- function(design) {
  design <- methods::as(design, "CsparseMatrix")
  vapply(seq_len(ncol(design)), function(i) {
    paste(design@i[seq.int(design@p[i] + 1, length.out = design@p[i + 1] -
                             design@p[i])], collapse = ",")
  }, character(1))
}
