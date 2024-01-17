source("./scripts/setup/setup.R")

ns <- c(50, 100, 200)
quantiles <- "zs"
method <- "quantile"

plot_res <- list()
for (j in 1:length(ns)) {

  n <- ns[j]

  sparse_beta <- c(rep(-2, 5), rep(-1, 5), rep(2, 5), rep(1, 5), rep(0, 80))
  dat <- gen_data(n = n, p = 100, beta = sparse_beta)

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
    boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, quantiles = "disturbed")
    res[[i]] <- ci.boot.ncvreg(boot_res, method = method, original_data = dat) %>%
      dplyr::mutate(width = (upper - lower), lambda = lambda_seq[i]) %>%
      dplyr::select(variable, width, lambda, estimate, lower, upper)
    if (any(is.na(res[[i]]$lower))) stop("Something wrong")
  }

  truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))

  plot_data <- do.call("rbind", res) %>%
    left_join(truth)

  plot_res[[j]] <- list("plot_data" = plot_data, "lambda_min" = lambda_min, "n" = n)

}

save(plot_res, file = glue("./rds/across_lambda_coverage_sparse_{quantiles}_{method}.rds"))

