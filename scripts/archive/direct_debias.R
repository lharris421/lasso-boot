rm(list = ls())
devtools::load_all()
library(hdrm)
library(ncvreg)

p <- 100
n <- 100
corr <- NULL
rho <- 0
SNR <- 1
laplace_beta <- rlaplace(p, rate = 1)
dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)
# dat <- gen_data_snr(n = 100, p = 100, beta = c(1, 2, 1, 0, 0, rep(0, 95)), rho = 0.5)
dat$y <- dat$y - mean(dat$y)
dat$X <- std(dat$X)
cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
lam <- cv_fit$lambda.min
lam_prop <- lam / max(cv_fit$lambda)

orig_betas <- coef(cv_fit$fit, lambda = cv_fit$lambda.min)
orig_int <- orig_betas[1]
orig_beta_est <- orig_betas[-1]

partial_resid <- (dat$y - orig_int) - (as.numeric(dat$X %*% orig_beta_est) - (dat$X * matrix(orig_beta_est, nrow=nrow(dat$X), ncol=ncol(dat$X), byrow=TRUE)))

var_num <- 2
zjs_orig <- (1/100)*colSums(dat$X * partial_resid)
bias1 <- (1/100) * t(dat$X[,var_num]) %*% dat$errs
bias2 <- (1/100) * t(dat$X[,var_num]) %*% dat$X[,-var_num] %*% (dat$beta[-var_num]*attr(dat$X, "scale")[-var_num] - orig_beta_est[-var_num])

(zjs_orig[var_num] - bias1 - bias2) / attr(dat$X, "scale")[var_num]
dat$beta[var_num]

corr_bias <- bias_est <- zjs <- betas <- lm_betas <- matrix(nrow = 10000, ncol = 100)
# bias_est <- numeric(1000)
# set.seed(1234)
for (i in 1:10000) {

  idx <- sample(1:100, replace = TRUE)
  newx <- dat$X[idx,]
  newx <- std(newx)
  newy <- dat$y[idx]
  newy <- newy - mean(newy)

  fit <- ncvreg(newx, newy, penalty = "lasso")
  curr_lambda <- lam_prop * max(fit$lambda)
  betas[i,] <- coef(fit, lambda = curr_lambda)[-1]

  (lm_betas[i,betas[i,] != 0] <- coef(lm(newy ~ -1 + newx[,betas[i,] != 0])))
  lm_betas[i,betas[i,] == 0] <- 0

  # resid <- newy - (newx %*% lm_betas[i,])

  partial_resid <- newy - (as.numeric(newx %*% betas[i,]) - (newx * matrix(betas[i,], nrow=nrow(newx), ncol=ncol(newx), byrow=TRUE)))

  zjs[i,] <- (1/100)*colSums(newx * partial_resid)

  tmp <- ifelse(betas[i,] == 0, 0, betas[i,] + (curr_lambda * sign(betas[i,]) - zjs[i,]))
  # bias_est[i] <- (1/100) * t(newx[,1,drop=FALSE]) %*% newx[,-1] %*% tmp[-1]

  for (j in 1:100) {
    bias_est[i,j] <- (1/100) * t(newx[,j,drop=FALSE]) %*% newx[,-j] %*% (lm_betas[i,-j] - betas[i,-j])
  }

  # bias_est[i] <- ifelse(betas[i,1] == 0, 0, curr_lambda * sign(betas[i,1]))

  if (i %% 100 == 0) {
    print(i)
    print(bias_est[i,1])
    # print(corr_bias[i,1])
  }

}
library(dplyr)


# samples <- zjs - bias_est - corr_bias
samples <- zjs - bias_est
# samples <- zjs
cis <- data.frame(t(apply(samples, 2, quantile, c(0.1, 0.9))))
names(cis) <- c("lower", "upper")
# cis <- cis %>%
#   mutate(bias_adj = colMeans(bias_est),
#          lower = lower - bias_adj, upper = upper - bias_adj)
mean(dat$beta >= cis$lower & dat$beta <= cis$upper)

# bias_est <- ifelse(betas == 0, 0, -lam * sign(betas))
# colMeans(bias_est)

# var <- 3
# (1/100) * t(dat$X[,var]) %*% dat$X[,-var] %*% colMeans(bias_est)[-var]

