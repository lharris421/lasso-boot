source("./scripts/setup.R")

dat <- gen_data_abn(n = 100, p = 100, a = 1, b = 1, rho = 0.5,
                    beta = c(2, rep(0, 99)))
boot_fit <- boot_ncvreg(dat$X, dat$y)

origX <- ncvreg::std(dat$X)

bias <- 2 -boot_fit$partial_correlations[,1]
#hist(bias)
abCorrs <- boot_fit$abCorrs

hist(abCorrs)
abline(v = (1/100)*origX[,1] %*% origX[,2], col = "red")

plot(abCorrs, bias)

summary(lm(bias ~ abCorrs))
