source("./scripts/setup/setup.R")

## Parameters
n <- 50
p <- 25
method <- "quantile"
methods <- c("ridge", methods)

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

cis <- list()
examples <- list()
for (i in 1:100) {

  print(i)
  dat <- gen_data_abn(n = n, p = p, a = 1, b = 1, rho = .99)

  for (j in 1:length(methods)) {
    if (i == 1) {
      cis[[j]] <- list()
    }
    if (methods[j] == "ridge") {
      ### Ridge
      ridge_cv <- cv_ridge(dat$X, dat$y)
      cis[[j]][[i]] <- cbind(confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min), method = "ridge")
    } else {
      ### Lasso-boot sample
      lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, max.iter = 1e9, quantiles = methods[j], nboot = nboot)
      cis[[j]][[i]] <- ci.boot.ncvreg(lassoboot, method = method) %>%
        dplyr::mutate(method = methods[j])
    }

    if (i == 37) {
      if (methods[j] == "ridge") {
        examples[[j]] <- cbind(confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))
      } else {
        examples[[j]] <- lassoboot
      }
    }

  }

}
names(cis) <- methods
names(examples) <- methods
save(cis, examples, file = glue("./rds/method_comparison_highcorr_100_{method}_n{n}_p{p}.rds"))
