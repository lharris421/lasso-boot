# Define the MCP penalty
mcp_penalty <- function(beta, lambda, gamma) {
  penalty <- ifelse(abs(beta) <= lambda * gamma,
                    lambda * abs(beta) - (beta^2 / (2 * gamma)),
                    (lambda^2 * gamma) / 2)
  sum(penalty)
}

# Define the log-posterior function
log_posterior <- function(beta, X, y, lambda, gamma) {
  # Likelihood (assume normal errors)
  residuals <- y - (X %*% beta)
  likelihood <- -0.5 * sum(residuals^2)

  # MCP penalty
  penalty <- -mcp_penalty(beta, lambda, gamma)

  # Log-posterior
  likelihood + penalty
}


# Metropolis-Hastings algorithm
metropolis_hastings <- function(X, y, lambda, gamma, iter = 1000000, burn_in = 100000, sigma_prop = 0.5) {

  accept <- numeric(iter - burn_in)
  n <- nrow(X)
  p <- ncol(X)

  # Initialize beta
  beta <- rep(0, p)

  # Store the samples
  samples <- matrix(NA, nrow = iter - burn_in, ncol = p)

  # Iterate
  for (i in 1:iter) {
    # Propose new beta from normal distribution
    beta_prop <- beta + rnorm(p, mean = 0, sd = sigma_prop)

    # Calculate acceptance ratio
    log_alpha <- log_posterior(beta_prop, X, y, lambda, gamma) - log_posterior(beta, X, y, lambda, gamma)
    alpha <- min(1, exp(log_alpha))

    # Accept or reject
    if (runif(1) < alpha) {
      beta <- beta_prop
      if (i > burn_in) accept[i - burn_in] <- 1
    }

    # Store the sample
    if (i > burn_in) samples[i - burn_in, ] <- beta
  }

  # Discard burn-in samples and return the rest
  return(list("samples" = samples, "acceptance_rate" = sum(accept / (iter - burn_in))))
}
