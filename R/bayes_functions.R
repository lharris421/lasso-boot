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

# ll <- function(beta, partial_residuals, sigma, xvar) {
#   sum(dnorm(partial_residuals - xvar * beta, mean = 0, sd = sigma, log = TRUE))
# }

density_function <- function(x, rate, partial_residuals, sigma, xvar) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
  return(exp(prior + llik))

}

# density_function <- function(x, rate, partial_residuals, sigma, xvar) {
#   prior <- log(dlaplace(x, rate = rate))
#   llik <- sapply(x, ll, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)
#   return(exp(prior + llik))
# }



obj <- function(beta, p, sigma, rate, xvar, partial_residuals, lower, upper) {

  denom <- integrate(
    density_function, lower = lower, upper = upper,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    subdivisions = 1000
  )$value

  if (p > .5) {
    prob <- integrate(
      density_function, lower = lower, upper = beta,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      subdivisions = 1000
    )$value * (1/denom)
  } else {
    prob <- integrate(
      density_function, lower = beta, upper = upper,
      rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
      subdivisions = 1000
    )$value * (1/denom)
    prob <- 1 - prob
  }

  return(p - prob)

}


post_quant <- function(p, sigma, rate, xvar, partial_residuals) {

  # Sample the function over the desired range
  ## tic(msg = "Getting range")
  beta_range <- seq(from = -10, to = 10, by = 0.01)  # Adjust these values as needed
  dens <- sapply(beta_range, density_function, rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar)

  # Determine the significant range
  non_zero_indices <- which(abs(dens) > 0)
  beta_min <- beta_range[min(non_zero_indices)]
  beta_max <- beta_range[max(non_zero_indices)]
  ## toc()

  ## tic(msg = "Determining Quantile")
  res <- uniroot(obj, c(beta_min, beta_max), p, sigma, rate, xvar, partial_residuals, beta_min, beta_max)
  ## toc()
  return(res$root)

}
