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
  tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))
  return(tmp)

}

density_function <- function(x, rate, partial_residuals, sigma, xvar) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik))

}

density_function_normalized <- function(x, rate, partial_residuals, sigma, xvar, normalizer) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik - log(normalizer)))

}

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

#
# post_quant <- function(p, post_mode, sigma, rate, xvar, partial_residuals) {
#
#
#   ## tic(msg = "Determining Quantile")
#   print(c(post_mode, sigma, rate))
#
#   print(".9")
#   print(obj2(post_mode, post_mode, .9, sigma, rate, xvar, partial_residuals))
#   # print(obj2(post_mode+30, post_mode, .9, sigma, rate, xvar, partial_residuals))
#   # print(".1")
#   # print(obj2(post_mode, post_mode, .1, sigma, rate, xvar, partial_residuals))
#   # print(obj2(post_mode-20, post_mode, .1, sigma, rate, xvar, partial_residuals))
#
#   if (p > .5) {
#     res <- uniroot(obj2, c(post_mode, post_mode + 30), post_mode, p, sigma, rate, xvar, partial_residuals)
#     print(p)
#     print(paste0("Results: ", res$root))
#   } else {
#     res <- uniroot(obj2, c(post_mode-20, post_mode), post_mode, p, sigma, rate, xvar, partial_residuals)
#     print(p)
#     print(paste0("Results: ", res$root))
#   }
#   ## toc()
#   return(res$root)
#
# }

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

  # xvals <- seq(-10, 10, by = .01)
  # yvals <- density_function_normalized(xvals, rate, partial_residuals, sigma, xvar, denom)
  # print(sum(yvals*.01))
  # plot(xvals, yvals, type = "l")
  # abline(v = 0, col = "red")

  res <- optim((1-p)/2, function(x) {
    q1 <- approx_quantile(x, rate, partial_residuals, sigma, xvar, denom)
    q2 <- approx_quantile(x + p, rate, partial_residuals, sigma, xvar, denom)
    q2 - q1
  }, method = "Brent", lower = 0, upper = 1 - p)
  c(approx_quantile(res$par, rate, partial_residuals, sigma, xvar, denom), approx_quantile(res$par + p, rate, partial_residuals, sigma, xvar, denom))
}




