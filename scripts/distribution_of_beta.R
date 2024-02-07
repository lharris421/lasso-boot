source("./scripts/setup/setup.R")

n <- 100
p <- 100
method <- "quantile"
all_res <- list()

## Function
accross_lambda_res <- function(dat, quantiles, method) {

  truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))

  if (length(unique(abs(truth$truth))) <= 4) {
    group_cutoffs <- sort(unique(abs(truth$truth)))[-1]
    group_names <- sort(unique(abs(truth$truth)))
  } else{
    group_cutoffs <- quantile(abs(truth$truth), c(.333, .667))
  }

  group_names <- c("Small", "Moderate", "Large")

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
    boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, quantiles = quantiles)
    res[[i]] <- ci.boot.ncvreg(boot_res, method = method) %>%
      dplyr::mutate(width = (upper - lower), lambda = lambda_seq[i]) %>%
      dplyr::select(variable, width, lambda, estimate, lower, upper) %>%
      dplyr::left_join(truth) %>%
      dplyr::mutate(
        covered = ifelse(truth >= lower & truth <= upper, 1, 0),
        group = group_names[(1 + (abs(truth) >= group_cutoffs[1]) + (abs(truth) >= group_cutoffs[2]))]
      )
    if (any(is.na(res[[i]]$lower))) stop("Something wrong")
  }

  return(list("plot_data" = do.call("rbind", res), "lambda_min" = lambda_min))

}

for (i in 1:length(methods)) {
  quantiles <- methods[i]

  ## Sparse
  set.seed(my_seed)
  sparse_beta <- c(rep(-2, 10), rep(-1, 10), rep(2, 10), rep(1, 10), rep(0, 60))
  dat <- gen_data(n = n, p = p, beta = sparse_beta)
  all_res[[1]] <- accross_lambda_res(dat, quantiles, method)

  ## Laplace
  rt <- 2
  set.seed(my_seed)
  dat <- gen_data(n = n, p = p, beta = rlaplace(p, rate = rt))
  all_res[[2]] <- accross_lambda_res(dat, quantiles, method)

  ## Normal
  set.seed(my_seed)
  dat <- gen_data(n = n, p = p, beta = rnorm(p, mean = 0, sd = 1))
  all_res[[3]] <- accross_lambda_res(dat, quantiles, method)

  ## T
  set.seed(my_seed)
  dat <- gen_data(n = n, p = p, beta = rt(p, 3))
  all_res[[4]] <- accross_lambda_res(dat, quantiles, method)

  save(all_res, file = glue("./rds/distribution_of_beta_n{n}_p{p}_{quantiles}.rds"))
}
