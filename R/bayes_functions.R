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
ll <- function(beta, partial_residuals, sigma, xvar) {

  # Compute the log-likelihood(s)
  tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))
  return(tmp)

}

## Computing density for a given beta for the posterior
density_function <- function(x, rate, partial_residuals, sigma, xvar, multiplier = 1) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik + sqrt(2)*length(partial_residuals) + log(multiplier)))

}

## Same as density above but takes an additional normalizing constant so that this density integrates to 1
density_function_normalized <- function(x, rate, partial_residuals, sigma, xvar, normalizer, multiplier = 1) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik + sqrt(2)*length(partial_residuals) - log(normalizer) + log(multiplier)))

}

## Objective for integrating from lower and upper bounds
obj <- function(beta, p, sigma, rate, xvar, partial_residuals, bounds, multiplier = 1) {


  ## Determine normalizing constant
  denom <- integrate(
    density_function, lower = bounds[1], upper = bounds[2],
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = multiplier
  )$value

  if (p > .5) {
    prob <- integrate(
      density_function_normalized, lower = bounds[1], upper = beta,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom, multiplier = multiplier
    )$value
  } else {
    prob <- integrate(
      density_function_normalized, lower = beta, upper = bounds[2],
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom, multiplier = multiplier
    )$value
    prob <- 1 - prob
  }

  return(p - prob)

}

post_quant <- function(p, post_mode, sigma, rate, xvar, partial_residuals) {

  ## One of the main hurdles is making sure this quantity is finite
  mltplyr <- 1 / density_function(post_mode, rate, partial_residuals, sigma, xvar, multiplier = 1)

  ## Determine the normalizing constant
  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  )$value

  ## Find the largest density, used for determining bounds of integration / plotting
  ymode <- density_function_normalized(post_mode, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

  ## Determine bounds
  step <- .01
  curr <- step
  while (TRUE) {
    xvals <- post_mode + c(-1, 1)*curr
    yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

    if (all(yvals < (ymode / 100))) {break} else {curr <- curr + step}
  }

  ## Redetermine normalizing constant based on new bounds
  denom <- integrate(
    density_function, lower = post_mode - curr, upper = post_mode + curr,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  )$value

  ## Determine the bounds
  res <- uniroot(obj, c(post_mode - curr, post_mode + curr), p, sigma, rate, xvar, partial_residuals, c(post_mode - curr, post_mode + curr))
  return(res$root)

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




