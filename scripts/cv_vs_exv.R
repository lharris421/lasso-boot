########################
## Setup ###############
########################
rm(list = ls())

library(progress) ## install.packages("progress")
library(glue)
library(ncvreg)

## Ridge CV Functions 
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

iterations <- 100
ridge_exv_sses <- lasso_exv_sses <- ridge_cv_sses <- lasso_cv_sses <- oracle_sses <- numeric(iterations)

## Data Scenarios
data_types <- list(
  list(p = 20, beta = c(2, -2, rep(0, 18))),
  list(p = 100, beta = c(1, -1, rep(0, 98))),
  list(p = 50, beta = c(rep(0.5, 8) * rep(c(-1, 1), each = 4), rep(0, 42))),
  list(p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60)))
)
data_type <- 1

pb <- progress_bar$new(
  format = "[:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 100
)
for (i in 1:iterations) {
  
  ##########
  ## Data ##
  ##########
  data <- hdrm::gen_data(n = 50, p = data_types[[data_type]]$p, beta = data_types[[data_type]]$beta)
  data_exv <- hdrm::gen_data(n = 50, p = data_types[[data_type]]$p, beta = data_types[[data_type]]$beta)
  
  ############
  ## Oracle ##
  ############
  oracle_beta <- coef(lm(data$y ~ -1 + data$X[,data$beta != 0]))
  
  ##############
  ## Using CV ##
  ##############
  lasso_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso")
  lasso_cv_beta <- coef(lasso_cv)[-1]
  ridge_cv <- cv_ridge(data$X, data$y)
  ridge_cv_beta <- ridge_cv$fit$beta[-1,ridge_cv$min]
  
  ###############
  ## Using EXV ##
  ###############
  lasso_exv_beta <- lasso_cv$fit$beta[-1,which.min(colMeans((data_exv$y - cbind(1, data_exv$X) %*% lasso_cv$fit$beta)^2))]
  ridge_exv_beta <- ridge_cv$fit$beta[-1,which.min(colMeans((data_exv$y - cbind(1, data_exv$X) %*% ridge_cv$fit$beta)^2))]
  
  
  ####################
  ## Calculate SSEs ##
  ####################
  oracle_sses[i] <- sum((data$beta[data$beta != 0] - oracle_beta)^2)
  lasso_cv_sses[i] <- sum((data$beta - lasso_cv_beta)^2)
  ridge_cv_sses[i] <- sum((data$beta - ridge_cv_beta)^2)
  lasso_exv_sses[i] <- sum((data$beta - lasso_exv_beta)^2)
  ridge_exv_sses[i] <- sum((data$beta - ridge_exv_beta)^2)
  pb$tick()
  
}

print(glue("Relative CV MSE of lasso to ridge: {round(mean(lasso_cv_sses) / mean(ridge_cv_sses), 3)}"))
print(glue("Relative Ex Validation MSE of lasso to ridge: {round(mean(lasso_exv_sses) / mean(ridge_exv_sses), 3)}"))
print(glue("Relative MSE of lasso CV to EXV: {round(mean(lasso_cv_sses) / mean(lasso_exv_sses), 3)}"))
print(glue("Relative MSE of ridge CV to EXV: {round(mean(ridge_cv_sses) / mean(ridge_exv_sses), 3)}"))
