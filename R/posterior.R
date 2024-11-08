#' Title
#'
#' @param X
#' @param y
#' @param cv_fit
#' @param lambda
#' @param sigma2
#' @param alpha
#' @param penalty
#' @param adjust_ss
#' @param studentize
#' @param relaxed
#'
#' @return
#' @export
#'
#' @examples
posterior <- function(X, y, cv_fit, lambda, sigma2, alpha = 0.05, penalty = "lasso",
                      adjust_ss = TRUE, relaxed = FALSE) {

  if (!missing(cv_fit) && class(cv_fit) != "cv.ncvreg") {
    stop("cv_fit must be an opbject of class cv.ncvreg.")
  }

  # Check if cv_fit$fit contains X and y, or if they are supplied
  if (missing(X) & missing(y) & missing(cv_fit)) {
    stop("You must supply X and y or cv_fit as an object of class cv.ncvreg")
  }

  if (missing(X)) {
    if (missing(cv_fit)) {
      stop("You must supply X or cv_fit as an object of class cv.ncvreg")
    } else if (missing(X) & !("X") %in% names(cv_fit$fit)) {
      stop("fit object in cv_fit missing X, please rerun cv.ncvreg with returnX = TRUE or supply X directly")
    }
  }

  if (missing(y) & missing(cv_fit)) {
    stop("You must supply y or cv_fit as an object of class cv.ncvreg")
  }

  if (!missing(cv_fit)) {
    if (!missing(y)) {
      message("Both cv_fit$fit and user-supplied y are provided; using cv_fit$fit's y.")
      y <- cv_fit$fit$y
    }
    if (!missing(X)) {
      message("Both cv_fit$fit and user-supplied X are provided; using cv_fit$fit's X.")
      X <- cv_fit$fit$X
    }
  }

  n <- nrow(X)
  X <- ncvreg::std(X)

  if (missing(lambda) & !missing(cv_fit)) {
    lambda <- cv_fit$lambda.min
  } else {
    cv_fit <- cv.ncvreg(X, y, penalty = penalty)
    lambda <- cv_fit$lambda.min
  }

  fit <- ncvreg(X, y, penalty = penalty)
  if (lambda < min(fit$lambda)) {
    fit <- ncvreg(X, y, penalty = penalty, lambda.min = 0.999*(lambda / max(fit$lambda)))
  } else if (lambda > max(fit$lambda)) {
    lambda <- max(fit$lambda)
  }
  bh_lambda <- coef(fit, lambda = lambda)
  rescale <- attr(X, "scale")
  bh_lambda <- bh_lambda[-1]

  intercept <- mean(y - as.numeric(X %*% bh_lambda))
  p <- length(bh_lambda)
  yhat <- intercept + as.numeric(X %*% bh_lambda)

  s <- bh_lambda != 0
  sh_lh <- sum(s)

  partial_residuals <- (y - yhat) + (X * matrix(bh_lambda, nrow = n, ncol = p, byrow=TRUE))
  b_bar <- (1/n)*colSums(X * partial_residuals)

  if (adjust_ss) {
    if (sh_lh > 0) {
      Xs <- X[,s,drop =FALSE]
      q_sh_default <- diag(n) - Xs %*% solve(t(Xs) %*% Xs, tol = 1e-12) %*% t(Xs)
    } else {
      q_sh_default <- diag(n)
    }

    ns <- numeric(p)
    for (j in 1:p) {

      if (!s[j]) {
        s_j <- s
        sh_j <- sh_lh
        q_sh <- q_sh_default
      } else {
        s_j <- s
        s_j[j] <- FALSE
        sh_j <- sum(s_j)

        if (sh_j > 0) {
          Xsj <- X[,s_j,drop=FALSE]
          q_sh <- diag(n) - Xsj %*% solve(t(Xsj) %*% Xsj, tol = 1e-12) %*% t(Xsj)
        } else {
          q_sh <- diag(n)
        }

      }

      ns[j] <- t(X[,j,drop=FALSE]) %*% q_sh %*% X[,j,drop=FALSE]


      if (relaxed) {
        b_bar[j] <- (t(X[,j,drop=FALSE]) %*% q_sh %*% y) / ns[j]
      }

    }
  }

  if (missing(sigma2)) sigma2 <- est_sigma(n = n, sh = sh_lh, resid = y - yhat)

  if (!adjust_ss) ns <- n
  ci <- ci_full_cond(b_bar, lambda = lambda, sigma2 = sigma2, n = ns,
                     alpha = alpha, penalty = penalty)

  data.frame(
    variable = colnames(X),
    estimate = bh_lambda / rescale,
    lower = ci$lower / rescale,
    upper = ci$upper / rescale,
    lambda = lambda,
    sigma = sqrt(sigma2)
  )

}
ci_full_cond <- function(z, lambda, sigma2, n, alpha, penalty = "lasso") {

  ## Tails being transferred on to (log probability in each tail)
  se <- sqrt(sigma2 / n)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)

  obs_p_lw <- obs_lw + ((z*lambda*n) / sigma2)
  obs_p_up <- obs_up - ((z*lambda*n) / sigma2)

  ## Find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))

  ps <- alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  lowers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )


  ps <- 1 - alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  uppers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )


  if (penalty == "MCP") {
    uppers <- sapply(uppers, firm_threshold_c, lambda, 3)
    lowers <- sapply(lowers, firm_threshold_c, lambda, 3)
  }
  return(data.frame(lower = lowers, upper = uppers, lambda = lambda))

}
soft_threshold <- function(z_j, lambda) {

  if (z_j > lambda) {
    return(z_j - lambda)
  } else if (abs(z_j) <= lambda) {
    return(0)
  } else if (z_j < -lambda) {
    return(z_j + lambda)
  }

}
## THIS IS TAKING A DRAW FROM THE POSTERIOR which is why we add back lambda
## z_j is not the best terminology here. Really z_j is beta j for this function
## should update to logically match
## probably should also update in case this was ever exactly 0???
firm_threshold_c <- function(z_j, lambda, gamma) {

  if (z_j != 0) {
    z_j <- z_j + sign(z_j)*lambda
  } else {
    stop("Logan you need to address this.")
  }

  if (abs(z_j) <= gamma*lambda) {
    return((gamma / (gamma - 1))*soft_threshold(z_j, lambda))
  } else {
    return(z_j)
  }

}

est_sigma <- function(resid, n, sh) {
  return((n - sh)^(-1) * sum(resid^2))
}
