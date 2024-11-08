## A script for testing the coverage rate for the posterior method (lasso / MCP penalties)
## Not really intended to save and veiw results, more just a general gut check

source("./scripts/setup.R")

p <- 100
n <- 50

true_rate <- sqrt(2*p)
true_lambda <- true_rate / n
penalty <- "lasso"

coverages <- numeric(100)
for (i in 1:100) {
  laplace_beta <- rlaplace(p, rate = 1)
  dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta)
  # dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = "autoregressive", rho = 0.8)
  
  # cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = penalty)
  
  coverages[i] <- 
    posterior(X = dat$X, y = dat$y, alpha = 0.2, adjust_ss = TRUE, relaxed = TRUE, lambda = true_lambda, sigma2 = 1) %>%
    # posterior(X = dat$X, y = dat$y, alpha = 0.2, adjust_ss = TRUE, lambda = true_lambda, sigma2 = 1) %>%
    #pipe_ncvreg(X = dat$X, y = dat$y, alpha = 0.2) %>%
    # posterior(X = dat$X, y = dat$y, lambda = true_lambda, sigma2 = 1, alpha = 0.2) %>%
    mutate(truth = dat$beta, covered = lower <= truth & truth <= upper) %>%
    pull(covered) %>%
    mean()
  
  print(coverages[i])
  
}
mean(coverages)


