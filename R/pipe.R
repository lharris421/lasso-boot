#' Title
#'
#' @param X
#' @param y
#' @param cv_fit
#' @param lambda
#' @param sigma
#' @param alpha
#' @param penalty
#' @param correction
#' @param original_n
#' @param inflate
#' @param relaxed
#'
#' @return
#' @export
#'
#' @examples
pipe_ncvreg <- function(
    X, y, cv_fit, lambda, sigma = NULL, alpha = 0.05, penalty = "lasso",
    original_n = FALSE, inflate = FALSE,
    relaxed = FALSE
) {

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

  if (missing(X)) {
    X <- cv_fit$fit$X
  }
  if (missing(y)) {
    y <- cv_fit$fit$y
  }

  ## Need to clean up logic, currently doesn't make sense
  if (missing(cv_fit)) {
    cv_fit <- cv.ncvreg(X, y, penalty = penalty)
    X <- cv_fit$fit$X
    y <- cv_fit$fit$y
  }

  if (missing(lambda)) {
    lambda <- cv_fit$lambda.min
  }

  n <- nrow(X)
  p <- ncol(X)

  bh_lambda <- coef(cv_fit$fit, lambda = lambda)
  rescale <- attr(X, "scale")
  bh_lambda <- bh_lambda[-1] * rescale
  intercept <- mean(y - as.numeric(X %*% bh_lambda))
  yhat <- intercept + as.numeric(X %*% bh_lambda)

  ## Same sigma for profile?? (should change name to conditional)
  s <- bh_lambda != 0
  sh_lh <- sum(s)
  if (is.null(sigma)) {
    sigma_h <- sqrt((n - sh_lh)^(-1) * sum((y - yhat)^2))
  } else {
    sigma_h <- sigma
  }

  partial_residuals <- (y - yhat) + (X * matrix(bh_lambda, nrow = n, ncol = p, byrow=TRUE))

  b_bar <- (1/n)*colSums(X * partial_residuals)

  if (sh_lh > 0) {
    Xs <- X[,s,drop =FALSE]
    p_sh_default <- Xs %*% solve(t(Xs) %*% Xs, tol = 1e-12) %*% t(Xs)
    q_sh_default <- diag(n) - p_sh_default
  } else {
    q_sh_default <- diag(n)
  }

  ses <- numeric(p)
  if (!original_n) {

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
      ses[j] <- sigma_h * sqrt(samp_size^-1)

      if (relaxed) {
        b_bar[j] <- (t(X[,j,drop=FALSE]) %*% q_sh %*% y) / samp_size
      }

    }

  } else {
    ses <- sigma_h / sqrt(n)
  }

  ts <- b_bar / ses
  ps <- 2 * (1 - pnorm(abs(ts)))
  qs <- p.adjust(ps, method = "BH")
  widths <- abs(qnorm(alpha / 2)) * ses

  data.frame(
    variable = colnames(X),
    estimate = b_bar / rescale,
    lower = (b_bar - widths) / rescale,
    upper = (b_bar + widths) / rescale,
    significance = qs,
    original_pvals = ps,
    lambda = lambda,
    sigma = sigma_h,
    betahat = bh_lambda / rescale,
    rescale_factor = rescale
  )

}
