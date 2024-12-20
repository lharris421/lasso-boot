#' Title
#'
#' @param X
#' @param y
#' @param fit
#' @param lambda
#' @param sigma
#' @param family
#' @param penalty
#' @param gamma
#' @param alpha
#' @param confidence_level
#' @param relaxed
#'
#' @return
#' @export
#'
#' @examples
pipe_bd <- function(X, y, fit, lambda, sigma,
                  family = c("gaussian", "binomial"),
                  penalty = c("MCP", "SCAD", "lasso"),
                  gamma = switch(penalty, SCAD = 3.7, 3),
                  alpha = 1, confidence_level = 0.95,
                  relaxed = FALSE
                  ) {

  ## If user wants more control, call ncvreg or cv.ncvreg directly
  family <- match.arg(family)
  penalty <- match.arg(penalty)

  ## I don't actually like this. Must supply X and y or fit object?
  if (missing(fit) & (missing(X) | missing(y))) {
    stop("Must supply a fit of class ncvreg or cv.ncvreg or both X and y.")
  }

  if ((!missing(fit) && is.null(fit$X)) & missing(X)) {
    stop("This procedure requires X. Either supply X, or fit the model using the option 'returnX = TRUE'")
  }

  ## Get a fit in some way, after this, every version of supplied X and y is the same
  if (missing(fit)) {
    if (missing(lambda)) {
      cv_fit <- cv.ncvreg(X, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, returnX = TRUE)
      fit <- cv_fit$fit
      lambda <- cv_fit$lambda.min
    } else {
      fit <- ncvreg(X, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, returnX = TRUE)
    }
  }

  if (is.null(fit$X)) {
    XX <- ncvreg::std(X)
  } else {
    XX <- fit$X
  }
  rescale_factorX <- attr(XX, "scale")
  yy <- fit$y ## not centered??, need to confirm always returned
  if (fit$family == "gaussian") yy <- yy - mean(yy) ## center y

  ## Check for nonsingular?? Can't do unless have original X
  # nonsingular <- attr(XX, "nonsingular")
  # if (length(nonsingular) != ncol(XX)) {
  #   warning("Columns in X are singular, inference will not be provided for these features")
  # }

  p <- ncol(XX)
  n <- nrow(XX)

  # load model (NEED TO UPDATE FOR BINOMIAL)
  if (fit$family == "gaussian") {
    beta <- coef(fit, lambda = lambda)[-1]
    beta <- beta * rescale_factorX
    S <- beta != 0
    S_hat <- which(S)
    N_hat <- which(!S)
  } else if (fit$family == "binomial") {
    beta <- coef(fit, lambda = lambda)
    S_hat <- which(beta[-1] != 0)
    N_hat <- which(beta[-1] == 0)
  }

  # initialize
  beta_PIPE <- numeric(p)
  sigma_PIPE <- numeric(p)

  if (fit$family == "gaussian"){

    # Compute sigma if not provided
    if (missing(sigma)) {
      sigma <- sqrt(crossprod(yy - XX %*% beta) / (n - sum(beta != 0)))
    }

    ## Compute pipe estimate
    yhat <- as.numeric(XX %*% beta)
    partial_residuals <- (yy - yhat) + (XX * matrix(beta, nrow = n, ncol = p, byrow = TRUE))
    if (!relaxed) beta_PIPE <- (1/n) * colSums(XX * partial_residuals)

    ## Compute PIPE Variance
    ## For features in the estimated support
    for (i in S_hat){

      S_hat_i <- S
      S_hat_i[i] <- FALSE

      ## Compute variance
      Xsi <- XX[,S_hat_i,drop = FALSE]
      Qsi <- diag(n) - tcrossprod(qr.Q(qr(Xsi)))
      adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qsi %*% XX[,i,drop = FALSE]
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)
      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop = FALSE]) %*% Qsi %*% yy) / adjusted_n
      }
    }

    # For null features
    if (length(S_hat) > 0 & length(N_hat) > 0) {
      Xs <- XX[,S_hat]
      Qs <- diag(n) - tcrossprod(qr.Q(qr(Xs)))
    } else {
      Qs <- diag(n)
    }

    for (i in N_hat) {

      adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qs %*% XX[,i,drop = FALSE]
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)

      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop = FALSE]) %*% Qs %*% yy) / adjusted_n
      }

    }

  } else if (fit$family == "binomial") {

    # compute pseudo outcome
    pii <- predict(fit, XX, type = "response", lambda = lambda)
    A <- pii * (1-pii)
    W <- diag(A) #
    v <- yy - pii # W and v change for different family of distribution
    Y_pse <- diag(1/A) %*% v + predict(fit, XX, type = "link", lambda = lambda)

    # compute weight matrix
    sqrtW <- sqrt(W)
    yy_pse <- sqrtW %*% Y_pse
    X_int <- cbind(rep(1,n), XX)

    # compute pipe for features in the estimated support
    for (i in S_hat) {

      S_hat_i <- beta[-1] != 0
      S_hat_i[i] <- FALSE

      Xs <- X_int[, c(TRUE, S_hat_i), drop = FALSE]
      XXs <- sqrtW %*% Xs
      xx <- crossprod(sqrtW, XX[,i])

      Qs <- diag(n) - tcrossprod(qr.Q(qr(XXs)))
      beta_PIPE[i] <- t(xx) %*% (yy - XXs %*% beta[c(TRUE, S_hat_i)]) / crossprod(xx)
      sigma_PIPE[i] <- sqrt(1 / (t(xx) %*% Qs %*% xx))

    }

    # compute pipe for null features
    S_hat <- beta[-1] != 0
    Xs <- X_int[,c(TRUE, S_hat), drop = FALSE]
    XXs <- sqrtW %*% Xs
    Qs <- diag(n) - tcrossprod(qr.Q(qr(XXs)))

    for (i in N_hat) {

      xx <- crossprod(sqrtW, XX[,i])
      beta_PIPE[i] <- t(xx) %*% (yy - XXs %*% beta[c(TRUE, S_hat)]) / crossprod(xx)
      sigma_PIPE[i] <- sqrt(1 / (t(xx) %*% Qs %*% xx))

    }

  }

  t <- beta_PIPE / sigma_PIPE
  pvalue <- (1 - pnorm(abs(t)))*2
  qvalue <- p.adjust(pvalue, method = "BH")
  ci_width <- qnorm(1 - ((1 - confidence_level)/2)) * sigma_PIPE

  if (fit$family == "binomial") beta <- beta[-1]

  res <- data.frame(
    variable = colnames(X),
    coef = beta / rescale_factorX,
    estimate = beta_PIPE / rescale_factorX,
    SE = sigma_PIPE,
    t = t,
    lower = (beta_PIPE - ci_width) / rescale_factorX,
    upper = (beta_PIPE + ci_width) / rescale_factorX,
    p.value = pvalue,
    p.adjust = qvalue,
    sigma = sigma,
    lamda = lambda
  )

  return(res)

}
