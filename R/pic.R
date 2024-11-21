pic <- function(
    X, y, penalty = "lasso", relaxed = TRUE
) {

  X <- ncvreg::std(X)
  n <- nrow(X)
  p <- ncol(X)

  cv_fit <- cv.ncvreg(X, y, penalty = penalty)

  logLiks <- numeric(length(cv_fit$lambda))
  eff_ss <- numeric(length(cv_fit$lambda))
  for (i in 1:length(cv_fit$lambda)) {

    bh_lambda <- coef(cv_fit$fit, lambda = cv_fit$lambda[i])[-1]
    yhat <- cv_fit$fit$linear.predictors[,i]

    ## Same sigma for profile?? (should change name to conditional)
    s <- bh_lambda != 0
    sh_lh <- sum(s)
    sigma_h <- sqrt((n - sh_lh)^(-1) * sum((y - yhat)^2))

    partial_residuals <- (y - yhat) + (X * matrix(bh_lambda, nrow = n, ncol = p, byrow=TRUE))
    b_bar <- (1/n)*colSums(X * partial_residuals)

    if (sh_lh > 0) {
      Xs <- X[,s,drop =FALSE]
      p_sh_default <- Xs %*% solve(t(Xs) %*% Xs, tol = 1e-12) %*% t(Xs)
      q_sh_default <- diag(n) - p_sh_default
    } else {
      q_sh_default <- diag(n)
    }

    logLik <- 0
    samp_sizes <- numeric(p)
    for (j in 1:p) {

      Xj <- X[,j,drop=FALSE]
      if (!s[j]) {

        s_j <- s
        sh_j <- sh_lh
        q_sh <- q_sh_default
        Xsj <- X[,s_j,drop=FALSE]

      } else {

        s_j <- s
        s_j[j] <- FALSE
        sh_j <- sum(s_j)

        if (sh_j > 0) {

          Xsj <- X[,s_j,drop=FALSE]
          Xsj2i <- solve(t(Xsj) %*% Xsj, tol = 1e-12)
          q_sh <- diag(n) - (Xsj %*% Xsj2i %*% t(Xsj))

        } else {

          Xsj <- X[,s_j,drop=FALSE]
          q_sh <- diag(n)

        }

      }

      ## May want to revisit
      samp_size <- t(X[,j,drop=FALSE]) %*% q_sh %*% X[,j,drop=FALSE]
      samp_sizes[j] <- samp_size
      if (relaxed) {
        b_bar[j] <- (t(X[,j,drop=FALSE]) %*% q_sh %*% y) / samp_size
      }
      logLik <- logLik + (samp_size / sigma_h^2) * (bh_lambda[j]^2 - 2*b_bar[j]*bh_lambda[j])


    }

    eff_ss[i] <- mean(samp_sizes)
    logLiks[i] <- logLik


  }

  return(logLiks + 2*colSums(cv_fit$fit$beta[-1,] != 0))

}
