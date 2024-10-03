rm(list = ls())
devtools::load_all()

library(MASS)  # For mvrnorm
library(coda)

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
traceplot(as.mcmc(samples))
