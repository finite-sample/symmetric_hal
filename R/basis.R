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

make_knots <- function(x, num_knots, knot_rule = "quantile") {
  knot_rule <- match.arg(knot_rule, c("quantile", "observed"))
  probs <- seq_len(num_knots) / (num_knots + 1)
  if (knot_rule == "observed") {
    ranks <- 1 + floor((nrow(x) - 1) * probs)
    ranks <- sort(unique(c(ranks, nrow(x) + 1 - ranks)))
    return(lapply(seq_len(ncol(x)), function(j) {
      unique(sort(x[, j])[ranks])
    }))
  }
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
  if (d == 1L) return(mains)
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

# Retain distinct functions even when their training columns agree.
make_dictionary <- function(x, num_knots = 3, symmetric = TRUE,
                            knot_rule = "quantile") {
  dictionary_features(make_knots(x, num_knots, knot_rule), symmetric)
}
