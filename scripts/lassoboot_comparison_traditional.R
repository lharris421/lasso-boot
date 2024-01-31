source("./scripts/setup/setup.R")

## Parameters
n <- 60
p <- a <- 60

## Data
beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 54, replace = TRUE))
dat <- gen_data(n = n, p = p, p1 = a, beta = beta, rho = 0)

cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
(lambda <- cv_fit$lambda.min)
(sigma2 <- cv_fit$cve[cv_fit$min])

res <- list()
for (i in 1:length(methods)) {
  res[[i]] <- boot.ncvreg(dat$X, dat$y, lambda = lambda, sigma2 = sigma2, verbose = TRUE, nboot = nboot, quantiles = methods[i])
}

save(dat, res, methods, file = glue("./rds/lassoboot_comparison_traditional_n{n}_p{p}.rds"))
