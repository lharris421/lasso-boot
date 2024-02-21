source("./scripts/setup/setup.R")

## Parameters
n <- 100
p <- a <- 100
SNR <- 1

## Data
res <- list()
current_seed <- my_seed
for (j in 1:100) {
  print(j)

  current_seed <- current_seed + j
  set.seed(current_seed)
  beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
  dat <- gen_data_snr(n = n, p = p, p1 = a, beta = beta, rho = 0, SNR = SNR)

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  print(lambda <- cv_fit$lambda.min)
  print(sigma2 <- cv_fit$cve[cv_fit$min])

  for (i in 1:length(methods)) {
    set.seed(current_seed)
    if (j == 1) {res[[i]] <- list()}
    tmp <- boot.ncvreg(dat$X, dat$y, lambda = lambda, sigma2 = sigma2, verbose = TRUE, nboot = nboot, method = methods[i], max.iter = 1e6)
    res[[i]][[j]] <- ci.boot.ncvreg(tmp, ci_method = ci_method) %>%
      mutate(method = methods[i], truth = dat$beta, group = j)
  }
}
names(res) <- methods

for (i in 1:length(methods)) {
  confidence_interval <- res[[methods[i]]]
  save(confidence_interval, file = glue("./rds/epsilon_conundrum_{methods[i]}_p{p}.rds"))
}

