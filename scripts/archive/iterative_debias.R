source("./scripts/setup/setup.R")

n <- p <- 100
corr <- NULL
rho <- 0
SNR <- 1

compute_partial_resid <- function(coefs, y, X) {
  (y - coefs[1]) - (
    as.numeric(X %*% coefs[-1]) - (X * matrix(coefs[-1], nrow = nrow(X), ncol = ncol(X), byrow=TRUE))
  )
}

compute_zs <- function(coefs, y, X) {
  (1/length(y)) * colSums(X * compute_partial_resid(coefs, y, X))
}

compute_mse <- function(coefs, y, X, truth) {
  mean((compute_zs(coefs, y, X) - truth)^2)
}

laplace_beta <- rlaplace(p, rate = 1)
dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)

cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")

boot_draws <- matrix(nrow = 1000, ncol = 100)
true_lambdas <- lambdas <- numeric(1000)

for (i in 1:1000) {

  idx <- sample(1:100, replace = TRUE)
  newy <- dat$y[idx]
  newy <- newy - mean(newy)
  newx <- ncvreg::std(dat$X[idx,])
  new_fit <- ncvreg(newx, newy, penalty = "lasso")
  lmin <- 0.999 * (cv_fit$lambda.min / max(new_fit$lambda))
  lmin <- ifelse(lmin < .05, lmin, .05)
  new_fit <- ncvreg(newx, newy, penalty = "lasso", lambda.min = lmin)
  tmp_coef <- coef(new_fit, lambda = cv_fit$lambda.min)
  true_lambdas[i] <- which.min(apply(coef(new_fit)[-1,] / attr(newx, "scale"), 2, function(x) mean((x - dat$beta)^2)))

  tmp_z <- compute_zs(tmp_coef, newy, newx)

  resids <- newy - (newx %*% tmp_z)
  newyy <- sample(resids, replace = TRUE) + (newx %*% tmp_z)
  newyy <- newyy - mean(newyy)
  tmp_fit <- ncvreg(newx, newyy, penalty = "lasso")

  which_lam <- which.min(apply(coef(tmp_fit), 2, compute_mse, as.numeric(newyy), newx, tmp_z))
  lambdas[i] <- which_lam
  which_lam <- true_lambdas[i]


  lmin <- 0.999* (tmp_fit$lambda[which_lam] / max(new_fit$lambda))
  lmin <- ifelse(lmin < .05, lmin, .05)
  new_fit <- ncvreg(newx, newy, penalty = "lasso", lambda.min = lmin)

  curr_lambda <- tmp_fit$lambda[which_lam]
  curr_lambda <- ifelse(curr_lambda > max(new_fit$lambda), max(new_fit$lambda), curr_lambda)
  updated_coefs <- coef(new_fit, lambda = curr_lambda)
  boot_draws[i,] <- compute_zs(updated_coefs, newy, newx) / attr(newx, "scale")
  if (i %% 100 == 0) print(i)

}


intervals <- data.frame(t(apply(boot_draws, 2, function(x) quantile(x, c(0.1, 0.9)))))
colnames(intervals) <- c("lower", "upper")

intervals %>%
  mutate(truth = dat$beta, covered = lower <= truth & truth <= upper) %>%
  pull(covered) %>%
  mean()

# Steps
# 1. Residual bootstrap dataset with z as beta estimates based on lambda cv
# 2. Fit model with data estimates
# 3. Calculate partial residuals and z_j est at every lambda
# 4. Select lambda with estimates relating to truth

