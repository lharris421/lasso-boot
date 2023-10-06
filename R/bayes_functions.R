################################################################################
## Laplace functions ###########################################################
################################################################################
dlaplace <- function(x, rate = 1) {
  dexp(abs(x), rate) / 2
}

qlaplace <- function(p, rate = 1) {
  if (p <= .5) {
    p <- (.5 - p)*2
    qexp(p, rate)*-1
  } else {
    p <- (p - .5)*2
    qexp(p, rate)
  }
}

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}
################################################################################

################################################################################
## Find Max Lambda #############################################################
################################################################################
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }
################################################################################

################################################################################
## Log-liklihood / posterior desnity function using partial residuals ##########
## (for Original method) #######################################################
################################################################################
ll <- function(beta, partial_residuals, sigma, xvar, ynew) {

  tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))

  return(tmp)

}

density_function <- function(x, rate, partial_residuals, sigma, xvar, normalizer = 1, multiplier = 0, log = FALSE, ynew) { # DONE

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, ynew = ynew)

  if (!log) {
    ret <- exp(prior + llik - log(normalizer) + multiplier)
  } else {
    ret <- prior + llik - log(normalizer) + multiplier
  }
  return(ret)

}
################################################################################

################################################################################
## Objective / function for posterior quantile finding #########################
################################################################################
obj <- function(beta, p, sigma, rate, xvar, partial_residuals, bounds, normalizer, multiplier, ynew) {

  prob <- integrate(
    density_function, lower = bounds[1], upper = beta,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    normalizer = normalizer, multiplier = multiplier, ynew = ynew
  )$value

  return(p - prob)

}

post_quant <- function(sig, post_mode, sigma, rate, xvar, partial_residuals, ynew) {

  ## Normalizing constant so that the log posterior = 0 at the posterior mode
  multiplier <- -density_function(
    x = post_mode, rate = rate, partial_residuals = partial_residuals,
    sigma = sigma, xvar = xvar, log = TRUE, ynew = ynew
  )

  ## Determine the normalizing constant (so now the posterior integrates to 1)
  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma,
    xvar = xvar, multiplier = multiplier, ynew = ynew
  )$value

  ## Find the largest density, used for determining bounds for uniroot
  ymode <- density_function(
    x = post_mode, rate = rate, partial_residuals = partial_residuals,
    sigma = sigma, xvar = xvar,
    normalizer = denom,
    multiplier = multiplier,
    ynew = ynew
  )

  ## Determine bounds
  step <- sigma
  curr <- step
  while (TRUE) {

    xvals <- post_mode + c(-1, 1)*curr
    yvals <- density_function(
      x = xvals, rate = rate, partial_residuals = partial_residuals,
      sigma = sigma, xvar = xvar,
      normalizer = denom,
      multiplier = multiplier,
      ynew = ynew
    )

    if (all(yvals < (ymode / 1000))) {
      break
    } else {
      curr <- curr + step
    }
  }

  ## redetermine
  denom <- integrate(
    density_function, lower = post_mode - curr, upper = post_mode + curr,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    multiplier = multiplier, ynew = ynew
  )$value

  ## Determine the bounds
  p_lower <- (1 - sig) / 2
  p_upper <- sig + p_lower

  lower <- uniroot(
    obj, c(post_mode - curr, post_mode + curr), p = p_lower, sigma = sigma,
    rate = rate, xvar = xvar, partial_residuals = partial_residuals,
    bounds = c(post_mode - curr, post_mode + curr), normalizer = denom,
    multiplier = multiplier, ynew = ynew
  )
  upper <- uniroot(
    obj, c(post_mode - curr, post_mode + curr), p = p_upper, sigma = sigma,
    rate = rate, xvar = xvar, partial_residuals = partial_residuals,
    bounds = c(post_mode - curr, post_mode + curr), normalizer = denom,
    multiplier = multiplier, ynew = ynew
  )
  return(c(lower$root, upper$root))

}
################################################################################

################################################################################
##### Bootstrapping / plotting / sim functions #################################
################################################################################
eb_boot <- function(beta, p = 60, b = 2, n = 100, nboot = 100, significance_level = .8, type = "univariate", lambda = "cv_once", prog = FALSE, sgm = 1, debias = FALSE, dat = NULL, time = FALSE, interval_type = "equal") {

  lower_p <- (1 - significance_level) / 2
  upper_p <- significance_level + lower_p

  if (time) tic(msg = "Overall")
  nlambda <- 100

  if (is.null(dat)) {

    if (time) tic(msg = "Generate data")
    dat <- genDataABN(beta = beta, p = p, a = length(beta), b = b, n = n, sgm = sgm)
    tbeta <- dat$beta

    X <- dat$X
    y <- dat$y
    p <- ncol(X)

    if (time) toc()

    if (lambda == "cv_once") {

      if (time) tic(msg = "Cross Validation")
      cv_res <- cv.ncvreg(X, y, penalty = "lasso")
      sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
      lam <- cv_res$lambda.min
      rate <- (lam*n / sigma2)
      if (time) toc()

    }


  } else {

    X <- dat$X
    y <- dat$y
    p <- ncol(X)
    n <- length(y)

    # if (lambda == "cv_once") {

      if (time) tic(msg = "Cross Validation")
      cv_res <- cv.ncvreg(X, y, penalty = "lasso")
      sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
      lam <- cv_res$lambda.min
      rate <- (lam*n / sigma2)
      if (time) toc()

      tbeta <- coef(cv_res$fit, lambda = lam)[-1]

    #}

  }

  lowers <- matrix(nrow = nboot, ncol = p)
  uppers <- matrix(nrow = nboot, ncol = p)

  if (prog) pb <- txtProgressBar(1, nboot, style = 3)

  if (time) tic(msg = "Bootstrapping")
  for (i in 1:nboot) {

    idx_new <- sample(1:length(y), replace = TRUE)
    ynew <- y[idx_new]
    xnew <- X[idx_new,,drop=FALSE]
    xnew <- ncvreg::std(xnew)

    if (lambda == "cv_once") {
      lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
      lambda_min <- lam - lam / 100 ## set min to be slightly smaller
      if (lambda_min > lambda_max | lam > lambda_max) {
        lambda_max <- lam + lam / 100
        nlambda <- 2
      }
      lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
      lasso_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    } else if (lambda == "cv_every") {
      cv_res <- cv.ncvreg(xnew, ynew, penalty = "lasso")
      sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
      lam <- cv_res$lambda.min
      coefs <- coef(cv_res$fit, lambda = lam)
      rate <- (lam*n / sigma2)
      lasso_fit <- ncvreg(xnew, ynew, penalty = "lasso")
    }

    coefs <- coef(lasso_fit, lambda = lam)

    if (type == "univariate") {

      ns_index <- attr(xnew, "nonsingular")
      post_modes <- coefs[-1]

      partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% post_modes) - (xnew * matrix(post_modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))

      z <- (1/n)*colSums(xnew * partial_residuals)
      se <- sqrt(sigma2 / n)

      obs_lw <- pnorm(0, z + lam, se)
      obs_up <- pnorm(0, z - lam, se, lower.tail = FALSE)

      enzls <- exp(n*z*lam / sigma2)^2

      lower_adj <- obs_lw + obs_up*enzls^(-1)
      upper_adj <- obs_up + obs_lw*enzls

      prop_lw <- obs_lw  / lower_adj

      lower <- ifelse(
        prop_lw >= lower_p,
        qnorm(lower_p * lower_adj, z + lam, se),
        qnorm(1 - ((1 - lower_p) * upper_adj), z - lam, se)
      )
      upper <- ifelse(
        prop_lw >= upper_p,
        qnorm(upper_p * lower_adj, z + lam, se),
        qnorm(1 - ((1 - upper_p) * upper_adj), z - lam, se)
      )

      rescale <- (attr(xnew, "scale")[ns_index])^(-1)
      lowers[i,ns_index] <- lower * rescale
      uppers[i,ns_index] <- upper * rescale

    } else {

      ## Beta specific
      for (j in attr(xnew, "nonsingular")) {

        idx <- as.numeric(which(attr(xnew, "nonsingular") == j))
        post_mode <- coefs[-1][idx]
        xvar <- xnew[,idx,drop=FALSE]

        partial_residuals <- ynew - (coefs[1] + xnew[,-idx,drop=FALSE] %*% coefs[-1][-idx])

        ## Need to compress this into a single call in future
        if (type == "original") {

          bounds <- (post_quant(significance_level, post_mode, sigma, rate, xvar, partial_residuals, ynew) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[idx])^(-1)

        } else if (type == "normal") {

          n <- length(ynew)
          norm_mean <- (2*sigma2*t(xvar) %*% partial_residuals) * (2*sigma2*t(xvar) %*% xvar + n^2*lam^2)^(-1)
          norm_var <- (2*sigma2^2) / (2*(t(xvar) %*% xvar)*sigma2 + n^2*lam^2)
          bounds <- (qnorm(c(lower_p, upper_p), norm_mean, sqrt(norm_var)) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[j])^(-1)

        } else if (type == "cadillac") {

          t2i_scale <- lam^2
          r <- ynew - (coefs[1] + xnew[,-idx] %*% coefs[-1][-idx])
          beta <- coefs[-1][idx]
          tau2i_mu <- sqrt((sigma2*lam^2) / beta^2)
          tau2 <- 1 / rinvgauss(1, tau2i_mu, t2i_scale)

          A <- (t(xnew[,idx, drop=FALSE]) %*% xnew[,idx,drop=FALSE]) + (1/tau2)
          mu <- solve(A)*(t(xnew[,idx,drop=FALSE]) %*% r)
          bounds <- (qnorm(c(lower_p, upper_p), mu, sqrt(sigma2*solve(A))) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[j])^(-1)

        } else {

          stop(paste0("Type: ", type, " not an option."))

        }

        uppers[i,j] <- bounds[2]
        lowers[i,j] <- bounds[1]

      }

    }

    if(prog) setTxtProgressBar(pb, i)

  }

  if (time) toc()
  if (time) toc() ## Overall

  return(list("lower" = lowers, "upper" = uppers, "truth" = tbeta))

}

plot_boot <- function(eb_boot, n = 30) {

  rm_lower <- apply(eb_boot[["lower"]], 2, function(x) sum(is.na(x))); names(rm_lower) <- names(eb_boot$truth)
  rm_upper <- apply(eb_boot[["upper"]], 2, function(x) sum(is.na(x))); names(rm_upper) <- names(eb_boot$truth)

  rm_lower <- rm_lower[rm_lower != 0]
  rm_upper <- rm_upper[rm_upper != 0]

  if (length(rm_lower) > 0 ) {
    print(paste0(length(rm_lower), " total variables with NA entries, summary: "))
    print(summary(rm_lower))
  }

  lowers <- apply(eb_boot[["lower"]], 2, mean, na.rm = TRUE)
  uppers <- apply(eb_boot[["upper"]], 2, mean, na.rm = TRUE)
  plot_res <- data.frame(truth = eb_boot[["truth"]], grp = names(eb_boot[["truth"]]), lower = lowers, upper = uppers) %>%
    dplyr::arrange(desc(abs(truth))) %>%
    head(n)
  plot_res$grp <- factor(plot_res$grp, levels = rev(plot_res$grp))

  print(plot_res)

  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = grp)) +
    geom_point(aes(x = truth, y = grp)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")
}

boot_ci <- function(eb_boot) {

  rm_lower <- apply(eb_boot[["lower"]], 2, function(x) sum(is.na(x))); names(rm_lower) <- names(eb_boot$truth)
  rm_upper <- apply(eb_boot[["upper"]], 2, function(x) sum(is.na(x))); names(rm_upper) <- names(eb_boot$truth)

  rm_lower <- rm_lower[rm_lower != 0]
  rm_upper <- rm_upper[rm_upper != 0]

  lowers <- apply(eb_boot[["lower"]], 2, mean, na.rm = TRUE)
  uppers <- apply(eb_boot[["upper"]], 2, mean, na.rm = TRUE)

  plot_res <- data.frame(estimate = eb_boot[["truth"]], variable = names(eb_boot[["truth"]]), lower = lowers, upper = uppers, method = "Lasso Boot")

  return(plot_res)

}


eb_boot_sim <- function(beta, p = 60, b = 2, n = 100, nboot = 100, nsim = 100, type = "original", debias = FALSE, interval_type = "equal") {

  overall_cov <- numeric(nsim)
  indiv_cov <- matrix(nrow = nsim, ncol = p)

  pb <- txtProgressBar(1, nsim, style = 3)

  ## p, beta, n
  for (iter in 1:nsim) {

    res <- eb_boot(beta = beta, p = p, b = b, n = n, nboot = nboot, type = type, debias = debias, interval_type = interval_type)

    lower_int <- apply(res[["lower"]], 2, mean)
    upper_int <- apply(res[["upper"]], 2, mean)
    tbeta <- res[["truth"]]

    indiv_cov[iter,] <- tbeta >= lower_int & tbeta <= upper_int
    overall_cov[iter] <- mean(indiv_cov[iter,])

    setTxtProgressBar(pb, iter)

  }

  return(list("overall_cov" = overall_cov, "indiv_cov" = indiv_cov, "truth" = tbeta))

}
################################################################################
plot_ci_comparison <- function(ci_df, nvars = 30) {

  plot_vars <- list()
  for (current_method in unique(ci_df$method)) {
    plot_vars[[current_method]] <- ci_df %>%
      filter(method == current_method) %>%
      dplyr::arrange(desc(abs(estimate))) %>%
      slice_head(n = nvars) %>%
      pull(variable)
  }
  plot_vars <- unique(unlist(plot_vars))

  plot_res <- ci_df %>%
    filter(variable %in% plot_vars) %>%
    dplyr::arrange(desc(abs(estimate)))

  plot_res$variable <- factor(plot_res$variable, levels = rev(plot_vars))

  gg <- plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable, color = method), alpha = .6) +
    geom_point(aes(x = estimate, y = variable, color = method), alpha = .6) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")

  return(gg)

}

