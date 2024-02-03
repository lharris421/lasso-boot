source("./scripts/setup/setup.R")

## Parameters
n <- 60
p <- a <- 60

## Data
res <- list()
current_seed <- my_seed
for (j in 1:100) {
  print(j)

  current_seed <- current_seed + j
  set.seed(current_seed)
  beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 54, replace = TRUE))
  dat <- gen_data(n = n, p = p, p1 = a, beta = beta, rho = 0)

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  (lambda <- cv_fit$lambda.min)
  (sigma2 <- cv_fit$cve[cv_fit$min])

  for (i in 1:length(methods)) {
    set.seed(current_seed)
    if (j == 1) {res[[i]] <- list()}
    tmp <- boot.ncvreg(dat$X, dat$y, lambda = lambda, sigma2 = sigma2, verbose = TRUE, nboot = nboot, quantiles = methods[i])
    res[[i]][[j]] <- ci.boot.ncvreg(tmp, method = method) %>%
      mutate(method = methods[i], truth = dat$beta, group = j)
  }
}
names(res) <- methods

for (i in 1:length(methods)) {
  confidence_interval <- res[[methods[i]]]
  save(confidence_interval, file = glue("./rds/epsilon_conundrum_{methods[i]}.rds"))
}

