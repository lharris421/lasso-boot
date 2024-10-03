devtools::load_all()
library(rstan)


stan_model_code <- "
data {
  int<lower=0> n;          // number of data points
  int<lower=0> p;          // number of predictors
  matrix[n, p] X;          // predictor matrix
  vector[n] y;             // response vector
  real<lower=0> lambda;    // penalty parameter
  real<lower=0> gamma;     // penalty parameter
}
parameters {
  vector[p] beta;          // regression coefficients
}
model {
  vector[n] mu;
  mu = X * beta;

  // Likelihood
  y ~ normal(mu, 1);   // Errors with standard deviation sigma

  // MCP penalty
  for (j in 1:p) {
    if (fabs(beta[j]) <= lambda * gamma) {
      target += -lambda * fabs(beta[j]) + 0.5 * beta[j]^2 / gamma;
    } else {
      target += -0.5 * lambda^2 * gamma;
    }
  }
}
"

# Example data
beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
dat <- gen_data_snr(n = 100, p = 100, p1 = 100, beta = beta, rho = 0, SNR = 1)
dat$X <- std(dat$X)
dat$y <- dat$y - mean(dat$y)
scalex <- drop(attr(dat$X, "scale")[1])
cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "MCP")
cv_fit$cve[cv_fit$min]

beta_est <- coef(cv_fit$fit, cv_fit$lambda.min)
int <- beta_est[1]
beta_est <- beta_est[-1]
partial_resid <- (dat$y - int) - (as.numeric(dat$X %*% beta_est) - (dat$X * matrix(beta_est, nrow=nrow(dat$X), ncol=ncol(dat$X), byrow=TRUE)))

# Data list for Stan
data_list <- list(n = nrow(dat$X), p = 1, X = dat$X[,1,drop=FALSE], y = partial_resid[,1], lambda = cv_fit$lambda.min / 100, gamma = 3 / 100) ## Marginal
# data_list <- list(n = nrow(dat$X), p = ncol(dat$X), X = dat$X, y = dat$y, lambda = cv_fit$lambda.min, gamma = 3) ## Full

# Compile and fit the model using Stan
fit <- stan(model_code = stan_model_code, data = data_list, iter = 100000, chains = 4, warmup = 90000, thin = 1)

# Check the results
print(fit)
plot(fit)
