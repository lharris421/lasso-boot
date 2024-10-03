source("./scripts/setup/setup.R")

p <- 100
n <- 100

coverages <- numeric(100)
for (i in 1:100) {
  laplace_beta <- rlaplace(p, rate = 1)
  dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta)

  coverages[i] <- confint(boot_ncvreg(dat$X, dat$y, verbose = FALSE)) %>%
    mutate(truth = dat$beta, covered = lower <= truth & truth <= upper) %>%
    pull(covered) %>%
    mean()

  print(coverages[i])

}
