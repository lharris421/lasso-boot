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

ll <- function(beta, partial_residuals, sigma, xvar) {

  # Compute the log-likelihood(s)
  ## (1/length(partial_residuals))
  tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))
  return(tmp)

}

density_function <- function(x, rate, partial_residuals, sigma, xvar, multiplier = 1) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik + sqrt(2)*length(partial_residuals) + log(multiplier)))

}

density_function_normalized <- function(x, rate, partial_residuals, sigma, xvar, normalizer, multiplier = 1) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik + sqrt(2)*length(partial_residuals) - log(normalizer) + log(multiplier)))

}

# Code / Objective for integrating from the mode
obj <- function(beta, p, sigma, rate, xvar, partial_residuals, bounds, multiplier = 1) {


  # denom <- integrate(
  #   density_function, lower = bounds[1], upper = bounds[2],
  #   rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = multiplier
  # )$value
  #
  # if (p > .5) {
  #   prob <- integrate(
  #     density_function_normalized, lower = bounds[1], upper = beta,
  #     rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
  #     normalizer = denom, multiplier = multiplier
  #   )$value
  # } else {
  #   prob <- integrate(
  #     density_function_normalized, lower = beta, upper = bounds[2],
  #     rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
  #     normalizer = denom, multiplier = multiplier
  #   )$value
  #   prob <- 1 - prob
  # }

  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = multiplier
  )$value

  if (p > .5) {
    prob <- integrate(
      density_function_normalized, lower = -Inf, upper = beta,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom, multiplier = multiplier
    )$value
  } else {
    prob <- integrate(
      density_function_normalized, lower = beta, upper = Inf,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      normalizer = denom, multiplier = multiplier
    )$value
    prob <- 1 - prob
  }



  return(p - prob)

}

# Code / Objective for integrating from the mode
# obj2 <- function(beta, post_mode, p, sigma, rate, xvar, partial_residuals) {
#
#
#   denom <- integrate(
#     density_function, lower = -Inf, upper = Inf,
#     rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar
#   )$value
#
#   xvals <- seq(-2, 2, by = .01)
#   yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom)
#   plot(xvals, yvals, type = "l")
#   abline(v = 0, col = "red")
#
#   if (p > .5) {
#     prob <- integrate(
#       density_function_normalized, lower = post_mode, upper = beta,
#       rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
#       normalizer = denom
#     )$value
#   } else {
#     prob <- integrate(
#       density_function_normalized, lower = beta, upper = post_mode,
#       rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
#       normalizer = denom
#     )$value
#   }
#
#   sig <- abs((1-2*p))
#   return((sig/2) - prob)
#
# }


post_quant <- function(p, post_mode, sigma, rate, xvar, partial_residuals, plot = FALSE) {

  ## One of the main hurdles
  mltplyr <- 1
  print(density_function(post_mode, rate, partial_residuals, sigma, xvar, multiplier = mltplyr))
  mltplyr <- 1 / density_function(post_mode, rate, partial_residuals, sigma, xvar, multiplier = mltplyr)
  print(density_function(post_mode, rate, partial_residuals, sigma, xvar, multiplier = mltplyr))

  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  )$value

  print(denom)


  ymode <- density_function_normalized(post_mode, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

  ## May want to revert back away from this if I can solve problem elsewhere
  step <- .01
  curr <- step
  while (TRUE) {
    xvals <- post_mode + c(-1, 1)*curr
    yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)

    if (all(yvals < (ymode / 100))) {break} else {curr <- curr + step}
  }

  # denom <- integrate(
  #   density_function, lower = post_mode - curr, upper = post_mode + curr,
  #   rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, multiplier = mltplyr
  # )$value
  #
  # print(denom)


  if (plot) {
    xvals <- seq(post_mode - curr, post_mode + curr, by = .01)
    yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom, multiplier = mltplyr)
    priory <- dlaplace(xvals, rate = rate) * (ymode / dlaplace(0, rate = rate))
    liky <- exp(ll(beta = xvals, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)) * length(partial_residuals)
    liky <- (liky / max(liky)) * ymode
    plot(xvals, yvals, type = "l")
    lines(xvals, priory, col = "blue", lty = 2)
    lines(xvals, liky, col = "blue", lty = 2)
    abline(v = post_mode, col = "red")
  }


  if (p > .5) {
    res <- uniroot(obj, c(post_mode - curr, post_mode + curr), p, sigma, rate, xvar, partial_residuals, c(post_mode - curr, post_mode + curr))
  } else {
    res <- uniroot(obj, c(post_mode - curr, post_mode + curr), p, sigma, rate, xvar, partial_residuals, c(post_mode - curr, post_mode + curr))
  }


  return(res$root)

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

# # Quantile function based on approximated CDF
# approx_quantile <- function(p, rate, partial_residuals, sigma, xvar, denom) {
#   uniroot(function(x) approx_cdf(x, rate, partial_residuals, sigma, xvar, denom) - p, interval = c(-10, 10))$root
# }
#
# # Optimize to find the narrowest interval containing p of the density
# find_narrowest_interval <- function(p, rate, partial_residuals, sigma, xvar) {
#
#   denom <- integrate(
#     density_function, lower = -10, upper = 10,
#     rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, subdivisions = 10000, rel.tol = 1e-12
#   )$value
#
#   # Plot
#   # xvals <- seq(-10, 10, by = .01)
#   # yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom)
#   # print(sum(yvals*.01)) ## Numerical integration
#   # plot(xvals, yvals, type = "l")
#   # abline(v = 0, col = "red")
#
#   res <- optim((1-p)/2, function(x) {
#     q1 <- approx_quantile(x, rate, partial_residuals, sigma, xvar, denom)
#     q2 <- approx_quantile(x + p, rate, partial_residuals, sigma, xvar, denom)
#     q2 - q1
#   }, method = "Brent", lower = 0, upper = 1 - p)
#   c(approx_quantile(res$par, rate, partial_residuals, sigma, xvar, denom), approx_quantile(res$par + p, rate, partial_residuals, sigma, xvar, denom))
# }




