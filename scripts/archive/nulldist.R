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

## Data arguments
nboot <- 1000
simulations <- 100

ideal_pvals <- pvals <- matrix(nrow = simulations, ncol = 100)

for (j in 1:simulations) {

  laplace_beta <- rlaplace(100, rate = 1)
  dat <- gen_data_snr(n = 100, p = 100, p1 = 100, beta = laplace_beta, SNR = 1)
  # dat <- hdrm::gen_data_abn()
  x_orig <- dat$X
  dat$X <- ncvreg::std(dat$X)
  rescale <- attr(dat$X, "scale")^(-1)

  dat$y <- dat$y - mean(dat$y)
  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min
  prop_lm <- lambda / max(cv_fit$lambda)

  orig_betas <- coef(cv_fit$fit, lambda = lambda)[-1]
  orig_zjs <- calc_zj(dat$X, dat$y, orig_betas)

  zjs <- matrix(nrow = 1000, ncol = ncol(dat$X))
  betas <- matrix(nrow = 1000, ncol = ncol(dat$X))

  for (i in 1:nboot) {

    idx <- sample(1:length(dat$y), replace = TRUE)
    xnew <- dat$X[idx,]
    xnew <- ncvreg::std(xnew)
    rescaleNew <- attr(xnew, "scale")^(-1)
    ynew <- dat$y[idx]
    ynew <- ynew - mean(ynew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda <- lambda_max * prop_lm
    lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
    if (lambda_min >= lambda_max | lambda >= lambda_max) { # Should review this
      lambda_max <- lambda + lambda / 100
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
    new_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    betanew <- coef(new_fit, lambda = lambda)[-1]

    zj <- calc_zj(xnew, ynew, betanew)
    zjs[i,] <- zj * rescaleNew

  }

  zjsnull <- matrix(nrow = 1000, ncol = ncol(dat$X))
  for (i in 1:nboot) {

    idx <- sample(1:length(dat$y), replace = FALSE)
    xnew <- dat$X
    xnew <- ncvreg::std(xnew)
    rescaleNew <- attr(xnew, "scale")^(-1)
    ynew <- dat$y[idx]
    ynew <- ynew - mean(ynew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda <- lambda_max * prop_lm
    lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
    if (lambda_min >= lambda_max | lambda >= lambda_max) { # Should review this
      lambda_max <- lambda + lambda / 100
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
    new_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    betanew <- coef(new_fit, lambda = lambda)[-1]

    zj <- calc_zj(xnew, ynew, betanew)
    zjsnull[i,] <- zj * rescaleNew


  }

  for (i in 1:ncol(zjs)) {

    pvals[j,i] <- mean(abs(zjsnull[,i]) > abs(orig_zjs[i]))

  }

  ideal_pvals[j,] <- 2*pnorm(abs(dat$beta) / sqrt(1/100), lower.tail = FALSE)

  print(j)
  print(mean(pvals[j,dat$beta == 0] < .2))
  print(mean(pvals[j,dat$beta != 0] < .2))

  plot(order(ideal_pvals[j,]), order(pvals[j,]))

}

