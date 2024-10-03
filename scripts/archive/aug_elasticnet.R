## Augmented attempt
source("./scripts/setup/setup.R")
dat <- hdrm::gen_data_abn()

x <- ncvreg::std(dat$X)
y <- dat$y - mean(dat$y)
p <- ncol(x)
n <- nrow(x)

alpha <- 0.1

cv_fit <- cv.glmnet(x, y, alpha = alpha, standardize = FALSE, thresh = 1e-18, maxit = 1e8)
lambda <- cv_fit$lambda.min

## Function to be used to find lambda max
find_thresh <- function(x, y, sample_size) { abs(t(x) %*% y) / sample_size }

## Univariate soft thresholding operator
univariate_soft_threshold <- function(z_j, lambda) {

  if (z_j > lambda) {
    return(z_j - lambda)
  } else if (abs(z_j) <= lambda) {
    return(0)
  } else if (z_j < -lambda) {
    return(z_j + lambda)
  } else {
    ## Me being lazy to handle errors
    print("Logan, you should take a closer look")
  }

}

## Check for convergence
convergence_check <- function(beta, beta_prev, tol = 1e-12) {

  all(abs((beta - beta_prev) / (beta_prev + tol)) < tol)

}

## Inner function only meant to be called within other functions
## Solves CD step for given lambda value when provided data, current lambda,
## and an initial value.
## This function assumes that X is standardized
lasso_cd_lambda <- function(X, y, curr_lambda, init, sample_size) {

  int <- mean(y)  ## Intercept w/ X standardized
  iter <- 0
  converged <- FALSE
  beta_prev <- beta <- init
  r_j <- y - (int + X %*% beta) ## Initial value of r_j
  n <- sample_size

  while (!converged & iter < 10000) {

    for (j in 1:ncol(X)) {

      x_curr <- X[,j]
      z_j <- c(((t(x_curr) %*% r_j)/n) + beta[j])
      beta[j] <- univariate_soft_threshold(z_j, curr_lambda)
      r_j <- r_j - (beta[j] - beta_prev[j]) * x_curr

    }

    iter <- iter + 1
    converged <- convergence_check(beta, beta_prev)
    beta_prev <- beta

  }

  return(beta)

}

lasso_cd <- function(X, y, sample_size, original_scale = TRUE) {

  lambda_max <- max(apply(X, 2, find_thresh, y, sample_size))
  lambda_min <- 1e-4 * lambda_max
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))

  ## Setup
  res <- matrix(nrow = ncol(X) + 1, ncol = length(lambda_seq))
  colnames(res) <- lambda_seq
  rownames(res) <- c("Intercept", colnames(X))

  for (i in 1:length(lambda_seq)) {

    if (i == 1) { init <- rep(0, ncol(X)) } else { init <- beta }

    curr_lambda <- lambda_seq[i]
    beta <- lasso_cd_lambda(X, y, curr_lambda, init, sample_size = sample_size)

    ## Save standardized or unstandardized estimates
    res[2:nrow(res),i] <- beta
    res[1,i] <- sum(y - (X %*% res[2:nrow(res),i])) / sample_size ## Re-solve for intercept

  }

  return(res)

}


coef_lasso <- function(res, lambdas, X, y, original_scale = TRUE) {

  original_lambdas <- as.numeric(colnames(res))
  tmp <- matrix(nrow = ncol(X) + 1, ncol = length(lambdas))
  colnames(tmp) <- lambdas
  rownames(tmp) <- c("Intercept", colnames(X))
  for (i in 1:length(lambdas)) {
    curr_lambda <- lambdas[i]
    if (curr_lambda >= max(lambdas)) {
      init <- rep(0, ncol(X))
    } else {
      ## Find nearest estimate already obtained upstream
      diffs <- original_lambdas - curr_lambda
      start <- tail(which(diffs >= 0), 1)
      init <- res[-1, start] * scales
    }
    beta <- lasso_cd_lambda(X, y, curr_lambda, init, sample_size = 100)

    tmp[2:nrow(tmp),i] <- beta
    tmp[1,i] <- mean(y - (X %*% tmp[2:nrow(tmp),i]))
  }
  return(tmp)
}


ynew <- c(y, rep(0, p))
xnew <- rbind(x, sqrt(n*(1 - alpha)*lambda)*diag(p))
xnew <- ncvreg::std(xnew)
rescale1 <- attr(xnew, "scale")
xnew <- xnew * sqrt(100 / 160)
# lambdanew <- max(apply(xnew, 2, find_thresh, ynew, 100)) * (lambda / max(apply(x, 2, find_thresh, y, 100)))
# final_adj <- (sd(xnew[,1]) / sd(x[,1]))^-1

tmp <- lasso_cd(xnew, ynew, sample_size = 100, original_scale = FALSE)

data.frame(
  # "aug" =  coef_lasso(tmp, alpha * lambda, xnew, ynew, original_scale = FALSE)[-1] * (sqrt(100 / 160) / rescale1),
  "aug" =  coef_lasso(tmp, alpha * lambda, xnew, ynew, original_scale = FALSE)[-1],
  "orig" = drop(coef(cv_fit$glmnet.fit, s = lambda, exact = TRUE, x = x, y = y))[-1]
) %>%
  mutate(diff = (aug - orig)^2) %>%
  pull(diff) %>%
  mean()

