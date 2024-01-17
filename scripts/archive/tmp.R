rm(list=ls())
#unloadNamespace("ncvreg")
#.libPaths("./local")
#library(ncvreg)
library(glmnet)
library(tictoc)
library(selectiveInference)
library(hdi)
library(ggplot2)
library(hdrm)
library(dplyr)
library(tidyr)
library(monomvn)


boot.ncvreg.2 <- function(X, y, cv_fit, lambda, sigma2, significance_level = 0.8, nboot = 100, ..., cluster, seed, returnCV=FALSE, verbose = TRUE, time = FALSE, quantiles = "sample") {

  if (time) tic(msg = "Overall")

  if (time) tic(msg = "Checks")
  if ((missing(X) | missing(y)) & (missing(cv_fit) || class(cv_fit) != "cv.ncvreg")) {
    stop("Either X and y or an object of class cv.ncvreg must be supplied.")
  }

  if (!missing(cv_fit)) {
    if (!all(c("X", "y") %in% names(cv_fit$fit))) {
      stop("fit object in cv_fit missing X and y, please rerun cv.ncvreg with returnX = TRUE or supply X and y directly (without also specifying a cv.ncvreg object)")
    }

    if(cv_fit$fit$penalty != "lasso") {
      stop(paste0("cv.ncvreg fit with ", cv_fit$fit$penalty, " penalty, but only 'lasso' penalty is currently supported for boot.ncvreg"))
    }

    if(cv_fit$fit$family != "gaussian") {
      stop(paste0("cv.ncvreg fit with ", cv_fit$fit$family, " family, but only 'gaussian' family is currently supported for boot.ncvreg"))
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
    rescale_original <- FALSE
  } else {
    rescale_original <- TRUE
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
    stop("Sorry, specification of alternative penality factors is not yet supported")
  }

  if (any(c("returnX", "warn", "convex") %in% names(args))) {
    warning(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("returnX", "warn", "convex")], collapse = ", "), " they are set to FALSE in cv.ncvreg"))
  }

  ## Note ignoring ncvreg arguments
  if ("family" %in% names(args)) {
    warning("Ignoring argument 'family', only guassian family is currently supported")
  }

  if ("penalty" %in% names(args)) {
    warning("Ignoring argument 'penalty', only lasso penalty is currently supported")
  }

  if (any(c("gamma", "alpha") %in% names(args))) {
    warning(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("gamma", "alpha")], collapse = " and "), ", not used for lasso penalty"))
  }

  if (time) toc()

  original_coefs <- NULL
  ## Will select lambda, won't estimate sigma^2 without selecting lambda
  if (missing(cv_fit)) {
    if (missing(lambda) | missing(sigma2)) {
      if (missing(lambda) & missing(sigma2)) {
        if (verbose) message("Using cross validation to select lambda and estimate variance")
      } else if (missing(lambda) & !missing(sigma2)) {
        if (missing(lambda) & verbose) message("Using cross validation to select lambda")
      } else if (!missing(lambda) & missing(sigma2) & verbose) {
        message("Using cross validation to estimate variance at supplied value of lambda using linear interpolation")
      }

      if (time) tic(msg = "Cross Validation")
      cv.args$penalty <- "lasso"
      cv.args$X <- X
      cv.args$y <- y
      if (!(missing(lambda))) {
        lambda_max <- max(apply(ncvreg::std(X), 2, find_thresh, y))
        lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
        nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
        if (lambda_min > lambda_max | lambda > lambda_max) {
          lambda_max <- lambda + lambda / 100
          nlambda <- 2
        }
        lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
        cv.args$lambda <- lambda_seq
      }
      if (!missing(cluster)) cv.args$cluster <- cluster ## NEED TO UPDATE
      cv_fit <- do.call("cv.ncvreg", c(cv.args, ncvreg.args))

      if (missing(lambda) & missing(sigma2)) {
        lambda <- cv_fit$lambda.min
        sigma2 <- cv_fit$cve[cv_fit$min]
      } else if (missing(lambda) & !missing(sigma2)) {
        lambda <- cv_fit$lambda.min
      } else if (!missing(lambda) & missing(sigma2)) {
        if (max(cv_fit$lambda) < lambda | min(cv_fit$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
        ## Make note about linear interpolation (or in documentation)
        ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
        l <- floor(ind)
        r <- ceiling(ind)
        w <- ind %% 1
        sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
      }

      original_coefs <- coef(cv_fit$fit, lambda = lambda)[-1]
      if (time) toc()
    }
  } else {
    if (time) tic(msg = "Setting arguments with supplied cv.ncvreg")
    if (!missing(sigma2) & verbose) message("Overriding variance estimate in cv.ncvreg object with user specified value for sigma2.")
    if (!missing(lambda) & verbose) message("Overriding selected lambda in cv.ncvreg object with user specified value for lambda.")
    if (!missing(lambda) & missing(sigma2)) warning("Estimating variance using CV but using user specified value of lambda. Estimated variance corresponds to interpolated CVE for supplied value of lambda.")
    if (missing(lambda)) lambda <- cv_fit$lambda.min
    if (missing(sigma2)) {
      if (max(cv_fit$lambda) < lambda | min(cv_fit$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
      ## Make note about linear interpolation (or in documentation)
      ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
    }
    X <- cv_fit$fit$X
    y <- cv_fit$fit$y
    original_coefs <- coef(cv_fit$fit, lambda = lambda)[-1]
    if (time) toc()
  }

  if (is.null(original_coefs)) {
    if (time) tic(msg = "Getting posterior mode")
    coef.args <- ncvreg.args
    coef.args$X <- X
    coef.args$y <- y
    coef.args$penalty <- "lasso"
    fit <- do.call("ncvreg", coef.args)
    original_coefs <- coef(fit, lambda = lambda)[-1]
    if (time) toc()
  }

  if (time) tic(msg = "Bootstrapping")
  modes <- matrix(nrow = nboot, ncol = ncol(X))
  per_draw <- 1
  draws <- matrix(nrow = nboot * per_draw, ncol = ncol(X))

  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "significance_level", "ncvreg.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(cl=cluster, X=1:nboot, fun=bootf.2, XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
  }

  for (i in 1:nboot) {
    if (!missing(cluster)) {
      res <- results[[i]]
    } else {
      res <- bootf.2(XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
    }
    draws[(1 + i*per_draw - per_draw):(i*per_draw),] <- res$draws
    modes[i,] <- res$modes
  }
  if (time) toc()

  val <- list(draws = draws, modes = modes, estimates = original_coefs, lamdba = lambda, sigma2 = sigma2)

  if (returnCV) val$cv.ncvreg <- cv_fit

  if (time) toc() ## For overall
  structure(val, class="boot.ncvreg.2")

}
bootf.2 <- function(XX, y, lambda, sigma2, significance_level = .8, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "sample") {

  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)

  modes <- numeric(p)
  if (time) toc()

  if (time) tic(msg = "Sample")
  idx_new <- sample(1:n, replace = TRUE)
  ynew <- y[idx_new]
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])
  nonsingular <- attr(xnew, "nonsingular")

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX

  if (time) toc()

  if (time) tic(msg = "Lambda Sequence")
  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }
  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  if (time) toc()

  if (time) tic(msg = "Fit ncvreg")
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq

  ## Ignores user specified lambda.min and nlambda
  # fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)

  coefs <- coef(fit, lambda = lambda)
  if (time) toc()

  if (time) tic(msg = "Compute Posterior")
  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))

  z <- (1/n)*colSums(xnew * partial_residuals)
  se <- sqrt(sigma2 / n)

  dmodes <- ifelse(
    modes <= 0,
    dnorm(modes, z + lambda, se, log = TRUE),
    dnorm(modes, z - lambda, se, log = TRUE)
  )


  spans <- runif(length(modes), 0, 3)

  accepted <- logical(length(modes))
  draws <- matrix(ncol = p, nrow = 1)
  iters <- 0
  while (any(!accepted) & iters < 1000) {
    for (i in 1:length(dmodes)) {
      if (!accepted[i]) {
        curr_sign <- sample(c(-1, 1), 1)
        curr_x <- modes[i] + curr_sign * spans[i] * se
        curr_thresh <- log(runif(1))


        curr_dens <- ifelse(
          curr_x <= 0,
          dnorm(curr_x, z[i] + lambda, se, log = TRUE) - dmodes[i],
          dnorm(curr_x, z[i] - lambda, se, log = TRUE) - dmodes[i]
        )

        if (curr_dens >= curr_thresh) {
          draws[1,i] <- curr_x
          accepted[i] <- TRUE
        } else if (sign(curr_x) != sign(modes[i])) {
          curr_x <- modes[i] + curr_sign * -1 * spans[i] * se
          curr_dens <- ifelse(
            curr_x <= 0,
            dnorm(curr_x, z[i] + lambda, se, log = TRUE) - dmodes[i],
            dnorm(curr_x, z[i] - lambda, se, log = TRUE) - dmodes[i]
          )
          if (curr_dens >= curr_thresh) {
            draws[1,i] <- curr_x
            accepted[i] <- TRUE
          }
        }

        spans[i] <- runif(1, 0, spans[i])
      }
    }
    iters <- iters + 1
  }

  draws <- draws * full_rescale_factor
  ## * full_rescale_factor

  if (time) toc()

  if (time) tic(msg = "Return result")

  modes[nonsingular] <- (modes * rescale) * rescaleX
  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA

  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)

}
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }


ci.boot.ncvreg.2 <- function(eb_boot, quiet = FALSE) {

  rm_draw <- apply(eb_boot[["draws"]], 2, function(x) sum(is.na(x))); names(rm_draw) <- names(eb_boot$estimates)
  rm_draw <- rm_draw[rm_draw != 0]

  all_draws <- eb_boot[["draws"]]
  lowers <- apply(all_draws, 2, function(x) quantile(x, 0.1, na.rm = TRUE))
  uppers <- apply(all_draws, 2, function(x) quantile(x, 0.9, na.rm = TRUE))

  ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")

  return(ci_info)

}


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

lasso_cis_s <- ridge_cis <- list()
for (i in 1:100) {

  dat <- gen_data_abn(n = 50, p = 25, a = 1, b = 1, rho = .99)

  ### Ridge
  ridge <- hdrm::ridge(dat$X, dat$y)
  ridge_cv <- cv_ridge(dat$X, dat$y)
  ridge_cis[[i]] <- confint(ridge, level = 0.8, lambda = ridge_cv$lambda.min)

  ### Lasso-boot sample
  lassoboot_s <- boot.ncvreg.2(dat$X, dat$y, verbose = FALSE, max.iter = 1e9, nboot = 1000)
  lasso_cis_s[[i]] <- ci.boot.ncvreg.2(lassoboot_s)

  if (i == 37) {
    lasso_example_s <- lassoboot_s
    ridge_example <- cbind(confint(ridge, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge, lambda = ridge_cv$lambda.min))
  }

}

ci.boot.ncvreg.2(lasso_example_s)

save(ridge_cis, lasso_cis_s,
     lasso_example_s, ridge_example, file = "./rds/method_comparison_highcorr_100_2.rds")
