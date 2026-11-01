# Rscript R/causal_path_check.R <completed_controls_dir>
invisible(loadNamespace("hal9001"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
input <- file.path(args[1], "hal9001_tight_original.rds")
saved <- readRDS(input)
stopifnot(is.null(saved$error), saved$method == "hal9001_tight")
data <- read.csv("data/lalonde_dw.csv")
data <- data[data$group %in% c("treated", "psid"), ]
x <- as.matrix(data[c("age", "educ", "black", "hisp", "married", "nodegree",
                      "re74", "re75")])
rows <- which(saved$folds$outer != 1 & data$treat == 1)
heldout <- which(saved$folds$outer == 1)
folds <- saved$folds$inner[[1]]$treated
diagnostics <- saved$diagnostics[["1_q1"]]$fit
selected <- match(diagnostics$lambda, diagnostics$lambda_grid)
stopifnot(!is.na(selected))
master <- hal9001::fit_hal(
  x[rows, ], data$re78[rows], family = "gaussian",
  lambda = diagnostics$lambda_grid[seq_len(selected)],
  fit_control = list(cv_select = FALSE, foldid = folds,
                     control = list(thresh = 1e-12, maxit = 1000000)),
  return_x_basis = TRUE
)
lambda <- tail(master$lasso_fit$lambda, 1)
stopifnot(master$lasso_fit$jerr == 0,
          abs(lambda / diagnostics$lambda - 1) < 1e-12)
design <- master$x_basis
beta <- as.vector(master$coefs[, ncol(master$coefs)])
y <- data$re78[rows]
residual <- beta[1] + as.vector(design %*% beta[-1]) - y
gradient <- as.vector(crossprod(design, residual)) / length(y)
kkt <- ifelse(beta[-1] != 0, abs(gradient + lambda * sign(beta[-1])),
              pmax(abs(gradient) - lambda, 0))
dense <- as.matrix(design)
centered <- sweep(dense, 2, colMeans(dense), `-`)
scale <- sqrt(colMeans(centered^2))
master_predictions <- predict(master, new_data = x[heldout, ])
master_predictions <- master_predictions[, ncol(master_predictions)]
train <- folds != 4
warnings <- character()
original <- withCallingHandlers(
  glmnet::glmnet(
    design[train, ], y[train], family = "gaussian", standardize = FALSE,
    lambda.min.ratio = 1e-4,
    control = list(thresh = 1e-12, maxit = 1000000)
  ),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
ratio <- 1e-4^(1 / 99)
last <- tail(original$lambda, 1)
stopifnot(original$jerr < 0, diagnostics$lambda < last)
extra <- ceiling(log(diagnostics$lambda / last) / log(ratio)) + 1
path <- c(original$lambda, last * ratio^seq_len(extra))
extended <- glmnet::glmnet(
  design[train, ], y[train], family = "gaussian", standardize = FALSE,
  lambda = path, control = list(thresh = 1e-12, maxit = 10000000)
)
stopifnot(extended$jerr == 0,
          tail(extended$lambda, 1) < diagnostics$lambda)
prediction <- function(fit, s) {
  as.vector(predict(fit, newx = design[!train, ], s = s, type = "response"))
}
old <- prediction(original, diagnostics$lambda)
endpoint <- prediction(original, last)
new <- prediction(extended, diagnostics$lambda)
stopifnot(identical(old, endpoint))
result <- data.frame(
  outer_fold = 1, component = "q1", inner_fold = 4,
  selected_lambda = diagnostics$lambda, master_jerr = master$lasso_fit$jerr,
  master_normalized_kkt = max(kkt[scale > 0] / scale[scale > 0]) / sd(y),
  master_scaled_intercept = abs(mean(residual)) / sd(y),
  master_prediction_movement = max(abs(
    master_predictions - saved$predictions$q1[heldout]
  )),
  original_inner_jerr = original$jerr, original_inner_last_lambda = last,
  extended_inner_jerr = extended$jerr,
  extended_inner_last_lambda = tail(extended$lambda, 1),
  clamp_movement = max(abs(old - endpoint)),
  extended_prediction_movement = max(abs(new - old)),
  original_validation_mse = mean((old - y[!train])^2),
  extended_validation_mse = mean((new - y[!train])^2)
)
readr::write_csv(result, "data/causal_path_check.csv")
hashes <- tools::md5sum(
  c("R/causal_path_check.R", "data/lalonde_dw.csv", input)
)
names(hashes) <- basename(names(hashes))
writeLines(trimws(c(capture.output(hashes), warnings,
                    capture.output(sessionInfo())), which = "right"),
           "data/causal_path_session.txt")
print(result, row.names = FALSE)
