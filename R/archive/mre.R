find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }

beta <- as.numeric(outer(c(0.5, 1, 2), c(-1, 1))); p <- 60; b <- 2; n <- 100
dat <- hdrm::genDataABN(beta = beta, p = p, a = length(beta), b = b, n = n)

X <- dat$X
y <- dat$y

cv_res <- cv.ncvreg(X, y, penalty = "lasso")
lam <- cv_res$lambda.min

for (i in 1:100) {

  idx_new <- sample(1:length(y), replace = TRUE)
  ynew <- y[idx_new]
  xnew <- X[idx_new,,drop=FALSE]
  xnew <- ncvreg::std(xnew)

  lambda_max <- max(apply(xnew, 2, find_thresh, ynew)); print(lambda_max)
  lambda_min <- lam - lam/100 ## set min to be slightly smaller
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))


  tryCatch({
    lasso_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda.min = lambda_min / lambda_max)
    coefs <- coef(lasso_fit, lambda = lam)
  },
  error = function(e) {
    print(e)
    print(lam)
    print(lasso_fit$lambda)
    }
  )

}

