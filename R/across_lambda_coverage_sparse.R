library(dplyr)
library(hdrm)

my_seed <- 189807771
ns <- c(50, 100, 150)
plot_res <- list()

for (j in 1:length(ns)) {

  n <- ns[j]

  rlaplace <- function(n, rate = 1) {
    rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
  }

  set.seed(my_seed)
  sparse_beta <- c(rep(-2, 5), rep(-1, 5), rep(2, 5), rep(1, 5), rep(0, 80))
  dat <- gen_data(n = n, p = 100, beta = sparse_beta)

  lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
  lambda_min <- lambda_max * 0.001
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = 0.001)
  lambda_min <- cv_fit$lambda.min

  if (!(lambda_max %in% lambda_seq)) {
    lambda_seq <- sort(c(lambda_seq, lambda_min), decreasing = TRUE)
  }

  res <- list()
  for (i in 1:length(lambda_seq)) {
    boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i])
    res[[i]] <- ncvreg:::ci.boot.ncvreg(boot_res) %>%
      dplyr::mutate(width = (upper - lower), lambda = lambda_seq[i]) %>%
      dplyr::select(variable, width, lambda, estimate, lower, upper)
  }

  truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))

  plot_data <- do.call("rbind", res) %>%
    left_join(truth)

  plot_res[[j]] <- list("plot_data" = plot_data, "lambda_min" = lambda_min, "n" = n)

}

save(plot_res, file = "./rds/across_lambda_coverage_sparse.rds")
