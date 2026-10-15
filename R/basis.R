# Zero-order HAL dictionaries of degree <= 2 on a shared knot grid.
#
# A feature is prod_{j in s} 1{o_j * (x_j - u_j) >= 0}: orientation +1 is
# x_j >= u_j, -1 is x_j <= u_j, so reflecting x_j -> -x_j maps the symmetric
# dictionary onto itself. The anchored dictionary uses +1 only (HAL's basis).
# The symmetric dictionary carries both orientations of main effects too.
# With a free intercept 1{x <= u} = 1 - 1{x > u}, which differs from
# 1{x >= u} only when u is a data value; with discrete covariates knots can
# be data values, and without the "-" main the reflected dictionary would
# lack the reflection of a "+" main.

make_knots <- function(x, num_knots) {
  probs <- seq_len(num_knots) / (num_knots + 1)
  lapply(seq_len(ncol(x)), function(j) {
    unique(unname(quantile(x[, j], probs, type = 7)))
  })
}

indicator <- function(xj, u, o) {
  if (o > 0) {
    outer(xj, u, `>=`) * 1
  } else {
    outer(xj, u, `<=`) * 1
  }
}

dictionary_features <- function(knots, symmetric) {
  d <- length(knots)
  main_orients <- if (symmetric) c(1L, -1L) else 1L
  mains <- do.call(rbind, lapply(main_orients, function(o) {
    do.call(rbind, lapply(seq_len(d), function(j) {
      data.frame(j = j, l = NA_integer_, u = knots[[j]], v = NA_real_,
                 oj = o, ol = NA_integer_)
    }))
  }))
  orients <- if (symmetric) {
    expand.grid(oj = c(1L, -1L), ol = c(1L, -1L))
  } else {
    data.frame(oj = 1L, ol = 1L)
  }
  pairs <- utils::combn(d, 2)
  inter <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(p) {
    j <- pairs[1, p]
    l <- pairs[2, p]
    grid <- expand.grid(v = knots[[l]], u = knots[[j]])
    do.call(rbind, lapply(seq_len(nrow(orients)), function(r) {
      data.frame(j = j, l = l, u = grid$u, v = grid$v,
                 oj = orients$oj[r], ol = orients$ol[r])
    }))
  }))
  rbind(mains, inter)
}

design_matrix <- function(x, features) {
  is_main <- is.na(features$l)
  cols <- vector("list", nrow(features))
  for (i in which(is_main)) {
    cols[[i]] <- indicator(x[, features$j[i]], features$u[i], features$oj[i])
  }
  inter <- which(!is_main)
  if (length(inter) > 0) {
    keys <- paste(features$j[inter], features$l[inter],
                  features$oj[inter], features$ol[inter])
    for (key in unique(keys)) {
      idx <- inter[keys == key]
      f <- features[idx[1], ]
      a <- indicator(x[, f$j], unique(features$u[idx]), f$oj)
      b <- indicator(x[, f$l], unique(features$v[idx]), f$ol)
      ua <- match(features$u[idx], unique(features$u[idx]))
      ub <- match(features$v[idx], unique(features$v[idx]))
      prod <- a[, ua, drop = FALSE] * b[, ub, drop = FALSE]
      for (m in seq_along(idx)) cols[[idx[m]]] <- prod[, m]
    }
  }
  Matrix::Matrix(do.call(cbind, cols), sparse = TRUE)
}

# Build the dictionary on training data. Interactions that are 1 on no more
# than 1/sqrt(n) of rows are dropped (hal9001's default screen, applied to
# every orientation alike); exact duplicates and constant columns are
# dropped. Main effects are never screened.
make_dictionary <- function(x, num_knots = 9, symmetric = FALSE) {
  features <- dictionary_features(make_knots(x, num_knots), symmetric)
  design <- design_matrix(x, features)
  share <- Matrix::colMeans(design)
  dup_of <- duplicated_to_first(design)
  keep <- dup_of == seq_along(dup_of) & share > 0 & share < 1 &
    (is.na(features$l) | share > 1 / sqrt(nrow(x)))
  kept <- features[keep, ]
  rownames(kept) <- NULL
  kept
}

duplicated_to_first <- function(design) {
  design <- methods::as(design, "CsparseMatrix")
  key <- vapply(seq_len(ncol(design)), function(i) {
    rows <- design@i[seq.int(design@p[i] + 1, length.out = design@p[i + 1] -
                               design@p[i])]
    paste(rows, collapse = ",")
  }, character(1))
  match(key, key)
}

is_interaction <- function(features) !is.na(features$l)

# hal9001 basis list for a dictionary, to be used with the design
# cbind(x, -x): a "-" indicator 1{x_j <= u} is hal9001's 1{-x_j >= -u} on
# column j + d. For the anchored dictionary only the first d columns are used.
to_hal9001 <- function(features, d) {
  col <- function(j, o) ifelse(o > 0, j, j + d)
  lapply(seq_len(nrow(features)), function(i) {
    f <- features[i, ]
    if (is.na(f$l)) {
      list(cols = col(f$j, f$oj), cutoffs = f$oj * f$u, orders = 0L)
    } else {
      list(cols = c(col(f$j, f$oj), col(f$l, f$ol)),
           cutoffs = c(f$oj * f$u, f$ol * f$v), orders = c(0L, 0L))
    }
  })
}
