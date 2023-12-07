library(monomvn) ## Other options EBglmnet, BGLR
## https://cs.gmu.edu/~pwang7/gibbsBLasso.html#:~:text=I%20implelemented%20a%20Gibbs%20sampler,priors%20for%20each%20regression%20coefficient.
library(hdrm)
library(ggplot2)
library(ncvreg)

# my_seed <- 189807771
# set.seed(my_seed)
dat <- genDataABN(n = 100, p = 10, a = 1, b = 1, rho = .99)

lassoboot <- boot.ncvreg.r(dat$X, dat$y, verbose = FALSE)
lambda2 <- lassoboot$lamdba^2
sigma2 <- lassoboot$sigma2
ci_lasso_boot <- ci.boot.ncvreg.r(lassoboot)

bayes_lasso <- blasso(dat$X, dat$y, T = 2000, lambda2 = lambda2, s2 = 1, RJ = FALSE, M = ncol(dat$X))
## Can set lambda as fixed with rd = FALSE

ci_bayes_lasso <- apply(bayes_lasso$beta[1001:2000,], 2, function(x) quantile(x, c(0.1, 0.5, 0.9)))


ci_bayes_lasso <- ci_bayes_lasso %>%
  t() %>%
  data.frame() %>%
  rename("lower" = X10.,"estimate" = X50., "upper" = X90.) %>%
  mutate(variable = ci_lasso_boot$variable, method = "Bayes Lasso") %>%
  select(estimate, variable, lower, upper, method)

all_methods <- list(ci_lasso_boot, ci_bayes_lasso)

do.call(rbind, all_methods) %>%
  ggplot() +
  geom_point(aes(x = estimate, y = variable, color = method)) +
  geom_errorbar(aes(xmin = lower, xmax = upper, y = variable, color = method)) +
  theme_bw() +
  labs(y = "Variable", x = "Estimate")

## cv.EBglmnet (only gives intervals for selected betas)
# library(EBglmnet)
# ebg <- cv.EBglmnet(dat$X, dat$y, prior = "lasso", nfolds = 10)
# ebg$fit

## gibbsBLasso (only provides estimates)
## gbl <- gibbsBLasso(dat$X, dat$y, max.steps = 100000)

## BGLR seems confusing
