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
calc_objective <- function(beta, X, y, lambda) {
  # obj <- (1/(2*length(y)))*crossprod(y - (X %*% beta)) - lambda*sum(abs(beta))
  # obj <- (1/(2*length(y)))*crossprod(y - (X %*% beta)) + lambda*sum(abs(beta))
  obj <- (1/(2*length(y)))*crossprod(y - (X %*% beta))
  return(obj)
}
## Data arguments
corr <- NULL
rho <- 0
p <- 100
n <- p * 1
nboot <- 1000
simulations <- 100
alpha <- .2
SNR <- 1

# laplace_beta <- rlaplace(p, rate = 1)
# dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)

coverages <- numeric(100)
for (j in 1:100) {

  dat <- hdrm::gen_data_abn()
  dat$X <- ncvreg::std(dat$X)
  rescale <- attr(dat$X, "scale")^(-1)

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min

  zjs <- matrix(nrow = 1000, ncol = ncol(dat$X))
  objs <- numeric(1000)

  for (i in 1:nboot) {

    idx <- sample(1:length(dat$y), replace = TRUE)
    xnew <- dat$X[idx,]
    xnew <- ncvreg::std(xnew)
    rescaleNew <- attr(xnew, "scale")^(-1)
    ynew <- dat$y[idx]
    ynew <- ynew - mean(ynew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
    if (lambda_min >= lambda_max | lambda >= lambda_max) { # Should review this
      lambda_max <- lambda + lambda / 100
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
    new_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    betanew <- coef(new_fit, lambda = lambda)[-1]

    zj <- calc_zj(xnew, ynew, betanew)

    adj_zj <- matrix(rnorm(1000 * ncol(xnew), 0, lambda), nrow = ncol(xnew)) + zj
    colsums <- colSums(abs(adj_zj))
    resample <- sample(1:length(colsums), replace = TRUE, prob = colsums)
    adj_zj <- t(adj_zj[,resample])

    obj_vals <- apply(adj_zj, 1, calc_objective, xnew, ynew, lambda)
    obj_vals <- 1 - (obj_vals / max(obj_vals))
    resample <- sample(1:length(obj_vals), replace = TRUE, prob = colsums)
    adj_zj <- adj_zj[resample,]

    zj <- colMeans(adj_zj)

    objs[i] <- drop(calc_objective(zj, xnew, ynew, lambda))
    zjs[i,] <- zj * (rescale * rescaleNew)

  }

  objs_adj <- 1 - (objs / max(objs))
  # objs_adj <- 1 - ((objs - min(objs)) / max(objs))
  # objs_adj <- objs
  resample <- sample(1:1000, replace = TRUE, prob = objs_adj)

  # ex <- 7
  #
  # hist(zjs[,ex])
  # hist(zjs[resample,ex])
  #
  # quantile(zjs[,ex], c(0.1, 0.5, 0.9))
  # quantile(zjs[resample,ex], c(0.1, 0.5, 0.9))

  # adjusted_samples <- zjs[resample,]
  adjusted_samples <- zjs
  lowers <- apply(adjusted_samples, 2, quantile, 0.05)
  uppers <- apply(adjusted_samples, 2, quantile, 0.95)

  res <- data.frame(truth = dat$beta, lower = lowers, upper = uppers)

  coverages[j] <- res %>%
    mutate(covered = lower <= truth & upper >= truth) %>%
    pull(covered) %>%
    mean()

  print(coverages[j])

}
