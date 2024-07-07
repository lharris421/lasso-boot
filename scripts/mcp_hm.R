library(MASS)  # For mvrnorm

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
metropolis_hastings <- function(X, y, lambda, gamma, iter = 1000000, burn_in = 100000) {

  accept <- numeric(iter)
  n <- nrow(X)
  p <- ncol(X)

  # Initialize beta
  beta <- rep(0, p)

  # Store the samples
  samples <- matrix(NA, nrow = iter, ncol = p)

  # Proposal standard deviation
  sigma_prop <- 0.5

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
      accept[i] <- 1
    }

    # Store the sample
    samples[i, ] <- beta
  }

  # Discard burn-in samples and return the rest
  return(list("samples" = samples[-(1:burn_in), ], "acceptance_rate" = sum(accept[-(1:burn_in)]) / (iter - burn_in)))
}

# Example usage
beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
dat <- gen_data_snr(n = 100, p = 100, p1 = 100, beta = beta, rho = 0, SNR = 1)
dat$X <- std(dat$X)
cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "MCP")

beta_est <- coef(cv_fit$fit, cv_fit$lambda.min)
int <- beta_est[1]
beta_est <- beta_est[-1]
partial_resid <- (dat$y - int) - (as.numeric(dat$X %*% beta_est) - (dat$X * matrix(beta_est, nrow=nrow(dat$X), ncol=ncol(dat$X), byrow=TRUE)))

samples <- metropolis_hastings(dat$X[,1,drop=FALSE], partial_resid[,1], cv_fit$lambda.min, gamma = 3)
print(samples$acceptance_rate)
samples <- samples$samples / attr(dat$X, "scale")[1]

# Check the samples
hist(samples)
quantile(samples, c(0.1, 0.9))
mean(quantile(samples, c(0.1, 0.9)))
library(coda)
traceplot(as.mcmc(samples))
