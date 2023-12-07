### ridge_comparison.R
### Bayes
bayes_lasso <- blasso(dat$X, dat$y, T = 2000, RJ = FALSE, M = ncol(dat$X), verb = 0)
ci_bayes_lasso <- apply(bayes_lasso$beta[1001:2000,], 2, function(x) quantile(x, c(0.1, 0.5, 0.9)))
ci_bayes_lasso <- ci_bayes_lasso %>%
  t() %>%
  data.frame() %>%
  rename("lower" = X10.,"estimate" = X50., "upper" = X90.) %>%
  mutate(variable = lasso_cis_s[[i]]$variable, method = "Bayes Lasso") %>%
  select(estimate, variable, lower, upper, method)
bayes_cis[[i]] <- ci_bayes_lasso

## method_comparison_lasso_sim.R
## Bayesian
start <- Sys.time()
bayes_lasso <- blasso(dat$X, dat$y, T = 2000, RJ = FALSE, M = ncol(dat$X), verb = 0)
ci_bayes_lasso <- apply(bayes_lasso$beta[1001:2000,], 2, function(x) quantile(x, c(0.1, 0.5, 0.9)))
ci_bayes_lasso <- ci_bayes_lasso %>%
  t() %>%
  data.frame() %>%
  rename("lower" = X10.,"estimate" = X50., "upper" = X90.) %>%
  mutate(variable = lassoboot_ci_s$variable, method = "Bayes Lasso") %>%
  select(estimate, variable, lower, upper, method)
bayes_lasso_coverage <- ci_bayes_lasso$lower <= laplace_beta & laplace_beta <= ci_bayes_lasso$upper
