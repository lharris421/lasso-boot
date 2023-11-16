library(glmnet)
library(tictoc)
library(selectiveInference)
library(hdi)
.libPaths("./local")
library(ncvreg)
.libPaths(.Library)
library(ggplot2)
library(hdrm)
library(dplyr)
library(tidyr)

my_seed <- 189807771
set.seed(my_seed)

cvf <- function(i, X, y, fold, cv.args) {
  XX <- X[fold!=i, , drop=FALSE]
  yy <- y[fold!=i]
  fit.i <- ridge(XX, yy)

  X2 <- X[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  loss <- ncvreg:::loss.ncvreg(y2, yhat, "gaussian")
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
cv_ridge <- function(X, y, nfolds = 10) {

    fit <- ridge(X, y)
    n <- length(y)
    E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))

    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds

    cv.args <- list()
    cv.args$lambda <- fit$lambda

    for (i in 1:nfolds) {
      res <- cvf(i, X, y, fold, cv.args)
      E[fold==i, 1:res$nl] <- res$loss
      Y[fold==i, 1:res$nl] <- res$yhat
    }

    ## Eliminate saturated lambda values, if any
    ind <- which(apply(is.finite(E), 2, all))
    E <- E[, ind, drop=FALSE]
    Y <- Y[, ind]
    lambda <- fit$lambda[ind]

    ## Return
    cve <- apply(E, 2, mean)
    cvse <- apply(E, 2, stats::sd) / sqrt(n)
    min <- which.min(cve)

    val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
                null.dev=mean(ncvreg:::loss.ncvreg(y, rep(mean(y), n), "gaussian")))
    return(val)

}


lasso_cis <- list()
ridge_cis <- list()
for (i in 1:100) {

  dat <- genDataABN(n = 100, p = 10, a = 1, b = 1, rho = .99)

  ### Ridge
  ridge <- hdrm::ridge(dat$X, dat$y)
  ridge_cv <- cv_ridge(dat$X, dat$y)
  ridge_cis[[i]] <- confint(ridge, level = 0.8, lambda = ridge_cv$lambda.min)

  ### Lasso-boot
  lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE)
  lasso_cis[[i]] <- ci.boot.ncvreg(lassoboot)

  if (i == 37) {
    lasso_example <- lassoboot
    ridge_example <- cbind(confint(ridge, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge, lambda = ridge_cv$lambda.min))
  }

}

save(ridge_cis, lasso_cis, lasso_example, ridge_example, file = "./rds/ridge_comparison.rds")
