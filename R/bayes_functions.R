## Laplace functions
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

## Log-liklihood using partial residuals
ll <- function(beta, partial_residuals, sigma, xvar, type = "full") {

  # Compute the log-likelihood(s)
  if (type == "univariate") {
    z <- (1/n)*as.numeric(t(xvar) %*% partial_residuals)
    tmp <- sapply(beta, function(x) dnorm(z - x, mean = 0, sd = sigma, log = TRUE))
  } else {
    tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))
  }

  return(tmp)

}

## Computing density for a given beta for the posterior
density_function <- function(x, rate, partial_residuals, sigma, xvar, normalizer = 1, multiplier = 1, type = "full", log = FALSE) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, type = type)

  ret <- ifelse(!log, exp(prior + llik - log(normalizer) + multiplier), prior + llik - log(normalizer) + multiplier)

  return(ret)

}

## Objective
obj <- function(beta, p, sigma, rate, xvar, partial_residuals, bounds, normalizer, multiplier = 1, type = "full") {

  prob <- integrate(
    density_function, lower = bounds[1], upper = beta,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    normalizer = normalizer, multiplier = multiplier, type = type
  )$value

  return(p - prob)

}

post_quant <- function(sig, post_mode, sigma, rate, xvar, partial_residuals) {

  ## Normalizing constant so that the log posterior = 0 at the posterior mode
  mltplyr <- -log_density_function(post_mode, rate, partial_residuals, sigma, xvar)

  ## Determine the normalizing constant (so now the posterior integrates to 1)
  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  )$value

  ## Find the largest density, used for determining bounds for uniroot
  ymode <- density_function_normalized(post_mode, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

  ## Determine bounds
  step <- sigma
  curr <- step
  while (TRUE) {
    xvals <- post_mode + c(-1, 1)*curr
    yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

    if (all(yvals < (ymode / 1000))) {
      break
    } else {
      curr <- curr + step
    }
  }

  ## redetermine
  denom <- integrate(
    density_function, lower = post_mode - curr, upper = post_mode + curr,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  )$value

  ## Determine the bounds
  p_lower <- (1 - sig) / 2
  p_upper <- sig + p_lower

  lower <- uniroot(obj, c(post_mode - curr, post_mode + curr), p_lower, sigma, rate, xvar, partial_residuals, c(post_mode - curr, post_mode + curr), denom, mltplyr)
  upper <- uniroot(obj, c(post_mode - curr, post_mode + curr), p_upper, sigma, rate, xvar, partial_residuals, c(post_mode - curr, post_mode + curr), denom, mltplyr)
  return(c(lower$root, upper$root))

}

plot_boot <- function(eb_boot) {

  lowers <- apply(eb_boot[["lower"]], 2, mean)
  uppers <- apply(eb_boot[["upper"]], 2, mean)
  plot_res <- data.frame(truth = eb_boot[["truth"]], grp = names(eb_boot[["truth"]]), lower = lowers, upper = uppers)

  plot_res %>%
    ggplot() +
    geom_point(aes(x = truth, y = grp)) +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = grp)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")

}

eb_boot <- function(beta, p = 60, b = 2, n = 100, niter = 100, plot = TRUE, type = "normal", prog = FALSE) {

  dat <- genDataABN(beta = beta, p = p, a = length(beta), b = b, n = n)

  tbeta <- dat$beta
  X <- dat$X
  y <- dat$y

  lowers <- matrix(nrow = niter, ncol = length(tbeta))
  uppers <- matrix(nrow = niter, ncol = length(tbeta))
  if (plot) {plot_data <- list()}

  if (prog) pb <- txtProgressBar(1, niter, style=3)

  for (i in 1:niter) {

    idx_new <- sample(1:length(y), replace = TRUE)
    ynew <- y[idx_new]
    xnew <- X[idx_new,,drop=FALSE]

    if (type == "univariate") {xnew <- ncvreg::std(xnew)}
    cv_res <- cv.ncvreg(xnew, ynew, penalty = "lasso")
    sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
    lam <- cv_res$lambda.min
    coefs <- coef(cv_res$fit, lambda = lam)
    rate <- (lam*n / sigma2)

    ## Beta specific
    for (j in 1:length(tbeta)) {

      post_mode <- coefs[-1][j]
      xvar <- xnew[,j,drop=FALSE]

      partial_residuals <- ynew - (coefs[1] + xnew[,-j,drop=FALSE] %*% coefs[-1][-j])

      ## Need to compress this into a single call in future
      if (type == "univariate") {
        z <- (1/n)*as.numeric(t(xvar) %*% partial_residuals)
        sigma <- sigma / sqrt(n)
        bounds <- post_quant_simp(.8, post_mode, sigma, rate, z) * attr(xnew, "scale")[j]^(-1)
        uppers[i,j] <- bounds[2]
        lowers[i,j] <- bounds[1]
      } else if (type == "normal") {
        n <- length(ynew)
        information_inv <- (2*sigma2^2) / (2*n*sigma2 + n^2*lam^2)
        score <- (1/sigma2)*(t(partial_residuals) %*% xvar - (t(xvar) %*% xvar + 1)*post_mode)
        norm_mean <- post_mode + (information_inv * score)
        uppers[i,j] <- qnorm(.9, norm_mean, sqrt(information_inv))
        lowers[i,j] <- qnorm(.1, norm_mean, sqrt(information_inv))
      } else if (type == "original") {
        bounds <- post_quant(.8, post_mode, sigma, rate, xvar, partial_residuals)
        uppers[i,j] <- bounds[2]
        lowers[i,j] <- bounds[1]
      } else if (type == "cadillac") {
        t2i_scale <- lam^2
        r <- ynew - (coefs[1] + xnew[,-j] %*% coefs[-1][-j])
        beta <- coefs[-1][j]
        tau2i_mu <- sqrt((sigma2*lam^2) / beta^2)
        tau2 <- 1 / rinvgauss(1, tau2i_mu, t2i_scale)

        A <- (t(xnew[,j, drop=FALSE]) %*% xnew[,j,drop=FALSE]) + (1/tau2)
        mu <- solve(A)*(t(xnew[,j,drop=FALSE]) %*% r)
        ci <- qnorm(c(.1, .9), mu, sqrt(sigma2*solve(A)))
        lowers[i,j] <- ci[1]; uppers[i,j] <- ci[2]
      } else {
        stop(paste0("Type: ", type, " not an option."))
      }

    }

    if(prog) setTxtProgressBar(pb, i)
  }

  return(list("lower" = lowers, "upper" = uppers, "truth" = tbeta))

}

eb_boot_sim <- function(beta, p = 60, b = 2, n = 100, niter = 100, nboot = 100, plot = FALSE) {

  overall_cov <- numeric(nboot)
  indiv_cov <- matrix(nrow = nboot, ncol = p)

  pb <- txtProgressBar(1, nboot, style=3)

  ## p, beta, n
  for (iter in 1:nboot) {

    res <- eb_boot(beta = beta, p = p, b = b, n = n, niter = niter, plot = plot)

    lower_int <- apply(res[["lower"]], 2, mean)
    upper_int <- apply(res[["upper"]], 2, mean)
    tbeta <- res[["truth"]]

    indiv_cov[iter,] <- tbeta >= lower_int & tbeta <= upper_int
    overall_cov[iter] <- mean(indiv_cov[iter,])

    setTxtProgressBar(pb, iter)

  }

  return(list("overall_cov" = overall_cov, "indiv_cov" = indiv_cov, "truth" = tbeta))

}

create_plots <- function(tmp) {

  plot_data <- tmp$plots

  # Initialize a list to store the plots
  plotly_plots <- list()

  # Iterate through the plot_data list
  for (i in seq_along(plot_data)) {
    # Extract the data for this plot
    data <- plot_data[[i]][[1]]

    # Create a data frame for ggplot
    plot_df <- data.frame(
      x = data$x,
      posterior = data$posterior,
      prior = data$prior,
      lik = data$lik
    )

    # Create the ggplot object
    p <- ggplot(plot_df, aes(x = x)) +
      geom_line(aes(y = posterior), color = "purple") +
      geom_line(aes(y = prior), color = "red") +
      geom_line(aes(y = lik), color = "blue") +
      geom_vline(xintercept = plot_data[[i]]$post_mode, color = "purple") +
      geom_vline(xintercept = plot_data[[i]]$lower, color = "purple", linetype="dashed") +
      geom_vline(xintercept = plot_data[[i]]$upper, color = "purple", linetype="dashed") +
      labs(title = paste("Lambda: ", round(plot_data[[i]]$lambda, 4)), x = "Beta", y = "Posterior Density")

    # Convert to ggplotly and store in the list
    plotly_plots[[i]] <- ggplotly(p)

  }

  return(plotly_plots)

}


################################################################################
######################### Currently Deprecated #################################
################################################################################

# Code / Objective for integrating from the mode
obj2 <- function(beta, post_mode, p, sigma, rate, xvar, partial_residuals) {


  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar
  )$value

  xvals <- seq(-2, 2, by = .01)
  yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom)
  plot(xvals, yvals, type = "l")
  abline(v = 0, col = "red")

  if (p > .5) {
    prob <- integrate(
      density_function_normalized, lower = post_mode, upper = beta,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom
    )$value
  } else {
    prob <- integrate(
      density_function_normalized, lower = beta, upper = post_mode,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom
    )$value
  }

  sig <- abs((1-2*p))
  return((sig/2) - prob)

}

post_quant2 <- function(p, post_mode, sigma, rate, xvar, partial_residuals) {


  if (p > .5) {
    res <- uniroot(obj2, c(post_mode, post_mode + 30), post_mode, p, sigma, rate, xvar, partial_residuals)
  } else {
    res <- uniroot(obj2, c(post_mode - 30, post_mode), post_mode, p, sigma, rate, xvar, partial_residuals)
  }
  return(res$root)

}



## Functions for finding HPD interval
# Approximate the CDF using numerical integration
approx_cdf <- function(x, rate, partial_residuals, sigma, xvar, denom) {
  res <- integrate(density_function_normalized, -10, x, rate, partial_residuals, sigma, xvar, denom, subdivisions = 10000, rel.tol = 1e-12)$value
  # print(paste0("Quantile: ", x))
  # print(res)
  return(res)
}

# Quantile function based on approximated CDF
approx_quantile <- function(p, rate, partial_residuals, sigma, xvar, denom) {
  uniroot(function(x) approx_cdf(x, rate, partial_residuals, sigma, xvar, denom) - p, interval = c(-10, 10))$root
}

# Optimize to find the narrowest interval containing p of the density
find_narrowest_interval <- function(p, rate, partial_residuals, sigma, xvar) {

  denom <- integrate(
    density_function, lower = -10, upper = 10,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, subdivisions = 10000, rel.tol = 1e-12
  )$value

  # Plot
  # xvals <- seq(-10, 10, by = .01)
  # yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom)
  # plot(xvals, yvals, type = "l")
  # abline(v = 0, col = "red")

  res <- optim((1-p)/2, function(x) {
    q1 <- approx_quantile(x, rate, partial_residuals, sigma, xvar, denom)
    q2 <- approx_quantile(x + p, rate, partial_residuals, sigma, xvar, denom)
    q2 - q1
  }, method = "Brent", lower = 0, upper = 1 - p)
  c(approx_quantile(res$par, rate, partial_residuals, sigma, xvar, denom), approx_quantile(res$par + p, rate, partial_residuals, sigma, xvar, denom))
}




