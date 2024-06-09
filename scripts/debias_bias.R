source("./scripts/setup/setup.R")
library(tictoc)

## Function
# Find partials
calc_zj <- function(X, y, beta) {
  partial_residuals <-  y - (as.numeric(X %*% beta) - (X * matrix(beta, nrow=nrow(X), ncol=ncol(X), byrow=TRUE)))
  z <- (1/length(y))*colSums(X * partial_residuals)
  return(z)
}
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }
biases <- matrix(nrow = 100, ncol = 100)
lambdas <- numeric(100)
cors <- numeric(100)
for (j in 1:100) {

  dat <- hdrm::gen_data_abn(p = 100, n = 100, a = 1, b = 2, rho = 0.5, SNR = 1)
  dat$X <- ncvreg::std(dat$X)
  rescale <- attr(dat$X, "scale")^(-1)

  dat$y <- dat$y - mean(dat$y)
  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min
  lambdas[j] <- lambda

  orig_betas <- coef(cv_fit$fit, lambda)[-1]
  orig_zjs <- calc_zj(dat$X, dat$y, orig_betas) * rescale

  biases[j,] <- orig_zjs - (dat$beta * rescale^(-1))

  if (sum(orig_betas != 0) > 1) {
    cors[j] <- pcor.test(dat$y, dat$X[,1], dat$X[,orig_betas!=0][,-1])$estimate
  } else {
    cors[j] <- 0
  }


}
