source("./scripts/setup.R")

dat <- gen_data_abn(b = 1)
cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
lambda <- cv_fit$lambda.min

scaledX <- ncvreg::std(dat$X)
orig_corr <- drop((1/100) * (scaledX[,1] %*% scaledX[,2]))

pb <- txtProgressBar(max = 1000)

corrs <- numeric(1000)
b2s <- numeric(1000)
for (i in 1:1000) {
  idx <- sample(1:100, replace = TRUE)
  scaledXnew <- ncvreg::std(dat$X[idx,])
  corrs[i] <- cor(scaledXnew[,1], scaledXnew[,2])
  
  b2s[i] <- coef(ncvreg(dat$X[idx,], dat$y[idx], penalty = "lasso"), lambda = lambda)[3]
  
  setTxtProgressBar(pb, i)
}

hist(corrs)
abline(v = orig_corr)
mean(corrs > orig_corr)

plot(b2s[b2s > 0], corrs[b2s > 0])
mean(corrs[b2s > 0])

