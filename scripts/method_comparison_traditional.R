source("./scripts/setup/setup.R")

## Parameters
n <- 60
nboot <- 1000
p <- a <- 60
quantiles <- "disturbed"

## Functions
# find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }
#
# trad_boot <- function(X, y, lambda) {
#
#   n <- length(y)
#   idx_new <- sample(1:n, replace = TRUE)
#   ynew <- y[idx_new]
#   xnew <- X[idx_new,,drop=FALSE]
#   xnew_std <- ncvreg::std(xnew)
#
#   lambda_max <- max(apply(xnew_std, 2, find_thresh, ynew))
#   lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
#   if (lambda_min > lambda_max | lambda > lambda_max) {
#     lambda_max <- lambda + lambda / 100
#     nlambda <- 2
#   }
#   lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
#   fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
#   res <- coef(fit, lambda = lambda)[-1]
#   return(res)
#
# }

## Data
beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 54, replace = TRUE))
dat <- gen_data(n = n, p = p, p1 = a, beta = beta, rho = 0)

cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
(lambda <- cv_fit$lambda.min)

trad_res <- boot.ncvreg(dat$X, dat$y, lambda = lambda, verbose = FALSE, nboot = nboot, quantiles = "mode")

## Ours
lasso_boot <- boot.ncvreg(dat$X, dat$y, lambda = lambda, verbose = FALSE, nboot = nboot, quantiles = quantiles)

save(dat, trad_res, lasso_boot, file = glue("./rds/method_comparison_traditional_{quantiles}_n{n}_p{p}.rds"))
