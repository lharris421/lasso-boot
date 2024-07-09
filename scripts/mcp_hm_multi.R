rm(list = ls())
library(MASS)  # For mvrnorm
library(tictoc)
library(coda)

# Example usage
beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
dat <- gen_data_snr(n = 100, p = 100, p1 = 100, beta = beta, rho = 0, SNR = 1)
dat$X <- std(dat$X)
cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "MCP")

tic()
samples <- metropolis_hastings(dat$X, dat$y, cv_fit$lambda.min, gamma = 3)
print(samples$acceptance_rate)
toc()

samples <- samples$samples

# Check the samples
scalex <- drop(attr(dat$X, "scale")[1])
hist(samples[,1] / scalex)
quantile(samples[,1] / scalex, c(0.1, 0.9))
mean(quantile(samples[,1] / scalex, c(0.1, 0.9)))

traceplot(as.mcmc(samples[,1]))
