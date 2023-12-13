rm(list=ls())
unloadNamespace("ncvreg")
.libPaths("./local")
library(ncvreg)
library(dplyr)
library(hdrm)

rt <- 2
ns <- c(30, 60, 90, 120)
plot_res <- list()

my_seed <- 189807771
set.seed(my_seed)

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}

for (j in 1:length(ns)) {

  n <- ns[j]

  true_lambda <- (1 / n) * rt

  laplace_beta <- rlaplace(60, rate = rt)
  dat <- gen_data(n = n, p = 60, beta = laplace_beta)

  lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
  lambda_min <- lambda_max * 0.001
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = 0.001)
  lambda_min <- cv_fit$lambda.min

  if (!any(abs(lambda_min - lambda_seq) < 1e-4)) {
    lambda_seq <- sort(c(lambda_seq, lambda_min), decreasing = TRUE)
  }

  res <- list()
  for (i in 1:length(lambda_seq)) {
    boot_res <- boot.ncvreg.r(X = dat$X, y = dat$y, lambda = lambda_seq[i])
    res[[i]] <- ci.boot.ncvreg.r(boot_res) %>%
      dplyr::mutate(width = (upper - lower), lambda = lambda_seq[i]) %>%
      dplyr::select(variable, width, lambda, estimate, lower, upper)
  }

  truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))

  plot_data <- do.call("rbind", res) %>%
    left_join(truth)

  plot_res[[j]] <- list("plot_data" = plot_data, "true_lambda" = true_lambda, "lambda_min" = lambda_min, "n" = n)

}

save(plot_res, file = "./rds/across_lambda_coverage_laplace.rds")
