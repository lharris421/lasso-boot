#' Bootstrap for ncvreg with gaussian outcomes
#'
#' Perform bootstrapping for the lasso at a single value of the
#' regularization parameter, lambda.
#'
#' At a specified value of lambda and variance
#' (which are selected / estimated by default if not provided)
#' the code does \code{nboot} iterations of taking a pairs bootstrap sample of
#' \code{X} and \code{y}, fitting \code{ncvreg}, and then uses this fit to
#' obtain hybrid bootstrap draws.
#'
#' Note that if a value of lambda is provided, but no value is provided for
#' sigma2 that sigma2 will be estimated using cross validation. Additionally,
#' the estimate for sigma2 will be interpolated if the supplied value of lambda
#' is not in the original sequence of lambda values.
#'
#' @param X The design matrix, without an intercept, as in \code{ncvreg}
#' @param y The response vector, as in \code{ncvreg}
#' @param cv_fit Instead of \code{X} and \code{y}, an object of type
#' \code{cv.ncvreg} with \code{returnX = TRUE} and \code{penalty = "lasso"}
#' @param penalty The penalty to be applied to the model.  Either "lasso" (the
#' default), "MCP", or "SCAD".
#' @param lambda A positive scalar value of type numeric, if left unspecified,
#' selected as \code{lambda.min} from \code{cv.ncvreg}
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details for \code{ncvreg}).
#' Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the
#' relative contributions from the lasso/MCP/SCAD penalty and the ridge, or L2
#' penalty.  \code{alpha=1} is equivalent to lasso/MCP/SCAD penalty, while
#' \code{alpha=0} would be equivalent to ridge regression.  However,
#' \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but
#' not exactly 0.
#' @param sigma2 A known / estimated value for the variance used in obtaining
#' bootstrap draws, if left unspecified, selected as the \code{cve}
#' corresponding \code{cv.ncvreg} to \code{lambda}
#' @param nboot The number of bootstrap iterations. Default is 1000.
#' @param ... Additional arguments to \code{ncvreg}/\code{cv.ncvreg}
#' @param cluster \code{cv.ncvreg} and the bootstrapping procedure can be run in
#' parallel across a cluster using the \code{parallel} package. The cluster must
#' be set up in advance using the \code{makeCluster} function from that package.
#' The cluster must then be passed to \code{boot_ncvreg} (see example).
#' @param seed You may set the seed of the random number generator in order to
#' obtain reproducible results. Note that the seed is set at the start of the
#' code an thus applies to \code{cv.ncvreg} as well as the bootstrap if
#' \code{cv.ncvreg} is run to select \code{lambda} and estimate \code{sigma2}.
#' If the user would like to specify a different seed to \code{cv.ncvreg}, they
#' should manually specify the call to \code{cv.ncvreg} in the arguments.
#' @param returnCV Whether to return \code{cv.ncvreg} object, \code{FALSE} by
#' default.
#' @param verbose Whether or not to print non-essential messages that highlight
#' potentially unanticipated behavior.
#' list(draws = draws, estimates = original_coefs, lambda = lambda, sigma2 = sigma2)
#' @return An object with S3 class \code{boot_ncvreg} containing: \describe {
#' \item{draws}{A \code{nboot} by \code{ncol(X)} matrix of the bootstrap draws}
#' \item{estimates}{A length \code{ncol(X)} vector of the estimates from the
#' lasso model fit on the original data corresponding to \code{lambda}.}
#' \item{cv.ncvreg}{If \code{returnCV == TRUE}, an object of type
#' \code{cv.ncvreg} that was fit to \code{X} and \code{y} to select
#' \code{lambda} and estimate \code{sigma2}. If user supplies \code{cv_fit},
#' will return the supplied \code{cv.ncvreg} object. If \code{lambda} and
#' \code{sigma2} are both specified, will return NULL.}}
#' @seealso \code{\link{ncvreg}}
#'
#' @examples
#'
#' data(Prostate)
#'
#' bootfit <- boot_ncvreg(Prostate$X, Prostate$y)
#'
#' ## Specifying cv_fit
#' bootfit <- boot_ncvreg(cv_fit = cv.ncvreg(Prostate$X, Prostate$y, penalty = "lasso", returnX = TRUE))
#'
#' ## Using a cluster
#' \dontrun{
#' library(parallel)
#' X <- Prostate$X
#' y <- Prostate$y
#' cl <- makeCluster(4)
#' bootfit <- boot_ncvreg(Prostate$X, Prostate$y, cluster = cl)}
#'
#' @export boot_ncvreg
boot_ncvreg <- function(X, y, cv_fit, penalty = "lasso",
                        lambda, gamma = switch(penalty, SCAD = 3.7, 3), alpha = 1,
                        sigma2_reed = FALSE, reselect_lambda = FALSE,
                        sigma2, nboot = 1000, ...,
                        cluster, seed, returnCV=FALSE, verbose = TRUE) {

  if ((missing(X) | missing(y)) & (missing(cv_fit) || class(cv_fit) != "cv.ncvreg")) {
    stop("Either X and y or an object of class cv.ncvreg must be supplied.")
  }

  if (!missing(cv_fit)) {
    if (!all(c("X", "y") %in% names(cv_fit$fit))) {
      stop("fit object in cv_fit missing X and y, please rerun cv.ncvreg with returnX = TRUE or supply X and y directly (without also specifying a cv.ncvreg object)")
    }

    if(cv_fit$fit$family != "gaussian") {
      stop(paste0("cv.ncvreg fit with ", cv_fit$fit$family, " family, but only 'gaussian' family is currently supported for boot_ncvreg"))
    }
  }

  if (!missing(y) & !missing(cv_fit)) {
    warning("Ignoring supplied values of y, using y from supplied cv.ncvreg object")
  }

  if (!missing(X) & !missing(cv_fit)) {
    if ("X" %in% names(cv_fit$fit)) {
      warning("Ignoring supplied values of X, using X from supplied cv.ncvreg object")
    }
  }

  if (!missing(X)) {
    rescale_original <- FALSE ## Supplied directly, will be scaled once in bootf
  } else {
    rescale_original <- TRUE ## Getting from cv.ncvreg where it has been scaled
  }

  # Coercion
  if (!missing(X)) {
    if (!is.matrix(X)) {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  }
  if (!missing(y)) {
    if (!is.double(y)) {
      tmp <- try(y <- as.double(y), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
    }
  }

  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }

  args <- list(...)
  if (length(args) > 0 & !missing(cv_fit)) {
    warning("Additional arguments are ignored when cv.ncvreg object supplied")
  }

  if (length(args) > 0 & any(is.null(names(args)))) {
    stop("Please supply names for all additional arguments passed to ...")
  }

  cv.args <- args[names(args) %in% c("nfolds", "fold", "returnY", "trace")]
  ncvreg.args <- args[names(args) %in% c("lambda.min", "nlambda", "eps", "max.iter", "dfmax")]
  if ("penalty.factor" %in% names(args)) {
    stop("Sorry, specification of alternative penality factors is not supported")
  }

  if (any(c("returnX", "warn", "convex") %in% names(args))) {
    warning(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("returnX", "warn", "convex")], collapse = ", "), " they are set to FALSE in cv.ncvreg"))
  }

  if ("family" %in% names(args)) {
    warning("Ignoring argument 'family', only guassian family is currently supported")
  }

  original_coefs <- NULL
  if (missing(cv_fit)) {
    if (missing(lambda) | missing(sigma2)) {

      if (missing(lambda) & missing(sigma2)) {
        if (verbose) message("Using cross validation to select lambda and estimate variance")
      } else if (missing(lambda) & !missing(sigma2)) {
        if (verbose) message("Using cross validation to select lambda")
      } else if (!missing(lambda) & missing(sigma2)) {
        if (verbose) message("Using cross validation to estimate variance at supplied value of lambda using linear interpolation")
      }

      cv.args$X <- X
      cv.args$y <- y
      cv.args$penalty <- penalty
      cv.args$alpha <- alpha
      cv.args$gamma <- gamma
      if (!(missing(lambda))) {

        nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
        lambda_seq <- ncvreg:::setupLambda(
          ncvreg::std(X), y, "gaussian", alpha = alpha,
          lambda.min = ifelse(nrow(X) > ncol(X), .001, .05),
          nlambda, penalty.factor=rep(1, ncol(X))
        )

        if (lambda < min(lambda_seq)) {

          lambda_min <- lambda - (lambda / 100)
          lambda_seq <- ncvreg:::setupLambda(
            ncvreg::std(X), y, "gaussian", alpha = alpha,
            lambda.min = lambda_min / max(lambda_seq),
            nlambda, penalty.factor=rep(1, ncol(X))
          )

        }

        if (lambda >= max(lambda_seq)) lambda <- max(lambda_seq) - (max(lambda_seq) / 100)

        cv.args$lambda <- lambda_seq

      }

      if (!missing(cluster)) cv.args$cluster <- cluster
      cv_fit <- do.call("cv.ncvreg", c(cv.args, ncvreg.args))

      if (missing(lambda) & missing(sigma2)) {
        lambda <- cv_fit$lambda.min
        sigma2 <- cv_fit$cve[cv_fit$min]
      } else if (missing(lambda) & !missing(sigma2)) {
        lambda <- cv_fit$lambda.min
      } else if (!missing(lambda) & missing(sigma2)) {
        ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
        l <- floor(ind)
        r <- ceiling(ind)
        w <- ind %% 1
        sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
      }

      fit <- cv_fit$fit
      which_lambda <- cv_fit$min
      original_coefs <- coef(fit, lambda = lambda)[-1]
    }
  } else {
    if (!missing(sigma2) & verbose) message("Overriding variance estimate in cv.ncvreg object with user specified value for sigma2.")
    if (!missing(lambda) & verbose) message("Overriding selected lambda in cv.ncvreg object with user specified value for lambda.")
    if (!missing(lambda) & missing(sigma2)) warning("Estimating variance using CV but using user specified value of lambda. Estimated variance corresponds to interpolated CVE for supplied value of lambda.")
    if (missing(lambda)) lambda <- cv_fit$lambda.min
    if (missing(sigma2)) {
      if (max(cv_fit$lambda) < lambda | min(cv_fit$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
      ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
    }
    fit <- cv_fit$fit
    X <- fit$X
    y <- fit$y
    which_lambda <- cv_fit$min
    original_coefs <- coef(fit, lambda = lambda)[-1]
  }

  if (is.null(original_coefs)) {

    coef.args <- ncvreg.args
    coef.args$X <- X
    coef.args$y <- y
    coef.args$penalty <- penalty
    coef.args$alpha <- alpha
    coef.args$gamma <- gamma

    nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
    lambda_seq <- ncvreg:::setupLambda(
      ncvreg::std(X), y, "gaussian", alpha = alpha,
      lambda.min = ifelse(nrow(X) > ncol(X), .001, .05),
      nlambda, penalty.factor = rep(1, ncol(X))
    )

    if (lambda < min(lambda_seq) | lambda > max(lambda_seq)) warning("Lambda outisde of range of original model fit, extending lambda sequence.")

    if (lambda < min(lambda_seq)) {

      lambda_min <- lambda - (lambda / 100)
      lambda_seq <- ncvreg:::setupLambda(
        ncvreg::std(X), y, "gaussian", alpha = alpha,
        lambda.min = lambda_min / max(lambda_seq),
        nlambda, penalty.factor=rep(1, ncol(X))
      )

    }

    if (lambda >= max(lambda_seq)) lambda <- max(lambda_seq) - (max(lambda_seq) / 100)

    coef.args$lambda <- lambda_seq

    fit <- do.call("ncvreg", coef.args)
    which_lambda <- which.min(abs(lambda_seq - lambda))
    original_coefs <- coef(fit, lambda = lambda)[-1]

  }

  partial_residuals <- y - (
    as.numeric(X %*% original_coefs) - (X * matrix(original_coefs, nrow = nrow(X), ncol = ncol(X), byrow=TRUE))
  )
  debiased_est <- numeric(ncol(X))
  nonsing <- attr(ncvreg::std(X), "nonsingular")
  debiased_est[nonsing] <- (1/nrow(X))*colSums(X * partial_residuals)


  debiased_draws <- bias <- point_estimates <- fc_draws <- partial_correlations <- matrix(nrow = nboot, ncol = ncol(X))


  if (sigma2_reed) {
    yhat <- fit$linear.predictors[,which_lambda]
    sh_lh <- sum(original_coefs != 0)
    sigma2 <- (length(y) - sh_lh)^(-1) * sum((y - yhat)^2)
  }


  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "ncvreg.args", "rescale_original", "penalty", "gamma", "alpha"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(
      cl=cluster, X=1:nboot, fun=bootf, XX = X, yy=y, lambda = lambda,
      sigma2 = sigma2, ncvreg.args = ncvreg.args, rescale_original = rescale_original,
      penalty = penalty, alpha = alpha, gamma = gamma
    )
  }

  for (i in 1:nboot) {
    if (!missing(cluster)) {
      res <- results[[i]]
    } else {
      res <- bootf(XX=X, yy=y, lambda = lambda, sigma2 = sigma2, reselect_lambda = reselect_lambda,
                   ncvreg.args = ncvreg.args, rescale_original = rescale_original,
                   penalty = penalty, alpha = alpha, gamma = gamma)
    }
    fc_draws[i,] <- res$fc_draws
    point_estimates[i,] <- res$point_estimates
    partial_correlations[i,] <- res$partial_correlations

  }

  colnames(fc_draws) <- names(original_coefs)
  colnames(point_estimates) <- names(original_coefs)
  colnames(partial_correlations) <- names(original_coefs)
  colnames(debiased_draws) <- names(original_coefs)
  names(debiased_est) <- names(original_coefs)
  val <- list(debiased_estimates = debiased_est,
              fc_draws = fc_draws, point_estimates = point_estimates,
              partial_correlations = partial_correlations,
              estimates = original_coefs, n = length(y),
              lambda = lambda, sigma2 = sigma2, penalty = penalty,
              alpha = gamma, gamma = gamma)

  if (returnCV) val$cv.ncvreg <- cv_fit

  structure(val, class="boot_ncvreg")

}
bootf <- function(XX, yy, lambda, sigma2, ncvreg.args, rescale_original = TRUE,
                  penalty = c("lasso", "MCP", "SCAD"), reselect_lambda = FALSE,
                  alpha = 1, gamma = switch(penalty, SCAD = 3.7, 3)) {

  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(yy)

  fc_draws <- point_estimates <- partial_correlations <- numeric(p)

  idx_new <- sample(1:n, replace = TRUE)
  ynew <- yy[idx_new]
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])

  nonsingular <- attr(xnew, "nonsingular")
  p_nonsingular <- length(nonsingular)

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleXX <- (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleXX <- 1
  }
  full_rescale_factor <- rescale * rescaleXX

  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- ncvreg:::setupLambda(
    ncvreg::std(xnew), ynew, "gaussian", alpha = alpha,
    nlambda, lambda.min = ifelse(n > p, .001, .05),
    penalty.factor=rep(1, ncol(xnew))
  )

  if (reselect_lambda) {
    cvnew <- cv.ncvreg(X = xnew, y = ynew, penalty = penalty, alpha = alpha, gamma = gamma)
    lambda <- cvnew$lambda.min
    sigma2 <- cvnew$cve[cvnew$min]
  }

  if (lambda < min(lambda_seq)) {
    lambda_min <- lambda - (lambda / 100)
    lambda_seq <- ncvreg:::setupLambda(
      ncvreg::std(xnew), ynew, "gaussian", alpha = alpha,
      lambda.min = lambda_min / max(lambda_seq), nlambda,
      penalty.factor=rep(1, ncol(xnew))
    )

  }

  if (lambda >= max(lambda_seq)) lambda <- max(lambda_seq) - (max(lambda_seq) / 100)

  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$lambda <- lambda_seq
  ncvreg.args$penalty <- penalty
  ncvreg.args$alpha <- alpha
  ncvreg.args$gamma <- gamma

  ## Ignore user specified lambda.min, nlambda since lambda sequence is being specified
  fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])

  modes <- coef(fit, lambda = lambda)[-1]

  if (alpha < 1) {
    ynew <- c(ynew, rep(0, p))
    xnew <- rbind(xnew, sqrt(n*(1 - alpha)*lambda)*diag(p))
    xnew <- ncvreg::std(xnew)
    xnew <- xnew * sqrt(n / (n + p))
    lambda <- lambda * alpha
  }

  partial_residuals <- ynew - (
    as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow = nrow(xnew), ncol = ncol(xnew), byrow=TRUE))
  )
  z <- (1/n)*colSums(xnew * partial_residuals)
  draws <- draw_full_cond(z, lambda, sigma2, n, p_nonsingular)

  if (penalty == "MCP") {
    draws <- sapply(draws, firm_threshold_c, lambda, gamma)
  } else if (penalty == "SCAD") {
    draws <- sapply(draws, scad_threshold_c, lambda, gamma)
  }

  fc_draws[nonsingular] <- draws * full_rescale_factor
  point_estimates[nonsingular] <- modes * full_rescale_factor
  partial_correlations[nonsingular] <- z * full_rescale_factor

  if (p_nonsingular < p) {
    fc_draws[!(1:p %in% nonsingular)] <- NA
    point_estimates[!(1:p %in% nonsingular)] <- NA
    partial_correlations[!(1:p %in% nonsingular)] <- NA
  }

  ret <- list(fc_draws, point_estimates, partial_correlations)
  names(ret) <- c("fc_draws", "point_estimates", "partial_correlations")
  return(ret)

}
draw_full_cond <- function(z, lambda, sigma2, n, p_nonsingular) {

  ## Tails being transferred on to (log probability in each tail)
  se <- sqrt(sigma2 / n)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)

  obs_p_lw <- obs_lw + ((z*lambda*n) / sigma2)
  obs_p_up <- obs_up - ((z*lambda*n) / sigma2)

  ## Find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))

  ps <- runif(p_nonsingular)

  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  draws <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )
  return(draws)

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
firm_threshold_c <- function(z_j, lambda, gamma) {

  z_j <- z_j + sign(z_j)*lambda

  if (abs(z_j) <= gamma*lambda) {
    return((gamma / (gamma - 1))*soft_threshold(z_j, lambda))
  } else {
    return(z_j)
  }

}
scad_threshold_c <- function(z_j, lambda, gamma) {

  z_j <- z_j + sign(z_j)*lambda

  if (abs(z_j) <= 2*lambda) {
    return(soft_threshold(z_j, lambda))
  } else if (abs(z_j) > 2*lambda & abs(z_j) <= gamma*lambda) {
    lambda_alt <- (gamma*lambda) / (gamma - 1)
    return(((gamma - 1) / (gamma - 2)) * soft_threshold(z_j, lambda_alt))
  } else {
    z_j
  }

}

