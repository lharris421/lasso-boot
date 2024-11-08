#' Title
#'
#' @param X
#' @param y
#' @param nfolds
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param X
#' @param y
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
ridge <- function(X, y, alpha) {
  ridge_cv <- cv_ridge(X, y)
  conf_ints <- confint(ridge_cv$fit, level = 1 - alpha, lambda = ridge_cv$lambda.min)
  cis <- bind_cols(conf_ints, "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))[-1,] %>%
    data.frame() %>%
    mutate(variable = colnames(X))
  colnames(cis) <- c("lower", "upper", "estimate", "variable")
  return(
    cis %>%
      select(variable, estimate, lower, upper) %>%
      mutate(
        lambda = ridge_cv$lambda.min,
        sigma = sqrt(ridge_cv$cve[ridge_cv$min])
      )
  )
}
