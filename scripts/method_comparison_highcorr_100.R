source("./scripts/setup/setup.R")

## Parameters
n <- 50
p <- 25
quantiles <- "disturbed"
method <- "quantile"

## Functions
cvf <- function(i, X, y, fold, cv.args) {
  XX <- X[fold!=i, , drop=FALSE]
  yy <- y[fold!=i]
  fit.i <- hdrm::ridge(XX, yy)

  X2 <- X[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  loss <- ncvreg:::loss.ncvreg(y2, yhat, "gaussian")
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
cv_ridge <- function(X, y, nfolds = 10) {

    ridge_fit <- hdrm::ridge(X, y)
    n <- length(y)
    E <- Y <- matrix(NA, nrow=n, ncol=length(ridge_fit$lambda))

    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds

    cv.args <- list()
    cv.args$lambda <- ridge_fit$lambda

    for (i in 1:nfolds) {
      res <- cvf(i, X, y, fold, cv.args)
      E[fold==i, 1:res$nl] <- res$loss
      Y[fold==i, 1:res$nl] <- res$yhat
    }

    ## Eliminate saturated lambda values, if any
    ind <- which(apply(is.finite(E), 2, all))
    E <- E[, ind, drop=FALSE]
    Y <- Y[, ind]
    lambda <- ridge_fit$lambda[ind]

    ## Return
    cve <- apply(E, 2, mean)
    cvse <- apply(E, 2, stats::sd) / sqrt(n)
    min <- which.min(cve)

    val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=ridge_fit, min=min, lambda.min=lambda[min],
                null.dev=mean(ncvreg:::loss.ncvreg(y, rep(mean(y), n), "gaussian")))
    return(val)

}

lasso_cis_s <- ridge_cis <- list()
for (i in 1:100) {

  dat <- gen_data_abn(n = n, p = p, a = 1, b = 1, rho = .99)

  ### Ridge
  ## ridge <- hdrm::ridge(dat$X, dat$y)
  ridge_cv <- cv_ridge(dat$X, dat$y)
  ridge_cis[[i]] <- confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min)

  ### Lasso-boot sample
  lassoboot_s <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, max.iter = 1e9, quantiles = quantiles, nboot = nboot)
  lasso_cis_s[[i]] <- ci.boot.ncvreg(lassoboot_s, method = method, original_data = dat)

  if (i == 37) {
    lasso_example_s <- lassoboot_s
    ridge_example <- cbind(confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))
  }

}

save(ridge_cis, lasso_cis_s,
     lasso_example_s, ridge_example, file = glue("./rds/method_comparison_highcorr_100_{quantiles}_{method}_n{n}_p{p}_span1.rds"))
