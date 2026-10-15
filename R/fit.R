# Each fitter returns list(predict = function(newx), beta, ...). `beta` is
# over the dictionary's features (NULL for tree and forest methods).

fit_cart <- function(x, y) {
  tree <- rpart::rpart(y ~ ., data = data.frame(y = y, x), method = "anova",
                       control = rpart::rpart.control(cp = 0, xval = 10))
  cp <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
  tree <- rpart::prune(tree, cp = cp)
  list(predict = function(newx) unname(predict(tree, data.frame(newx))),
       tree = tree)
}

fit_rf <- function(x, y) {
  forest <- ranger::ranger(y ~ ., data = data.frame(y = y, x),
                           num.threads = 1)
  list(predict = function(newx) {
    predict(forest, data.frame(newx), num.threads = 1)$predictions
  })
}

# HAL on a dictionary: lasso on unstandardized indicators, so the L1 norm
# of the coefficients is the dictionary penalty V of the fitted function.
fit_lasso <- function(design, y, features, relax = FALSE, nfolds = 10) {
  cv <- glmnet::cv.glmnet(design, y, alpha = 1, standardize = FALSE,
                          relax = relax, nfolds = nfolds)
  if (relax) {
    s <- cv$relaxed$lambda.min
    gamma <- cv$relaxed$gamma.min
  } else {
    s <- cv$lambda.min
    gamma <- 1
  }
  beta <- as.vector(coef(cv, s = s, gamma = gamma))[-1]
  list(
    predict = function(newx) {
      as.vector(predict(cv, design_matrix(newx, features), s = s,
                        gamma = gamma))
    },
    beta = beta
  )
}
