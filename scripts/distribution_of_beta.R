source("./scripts/setup/setup.R")

n <- 100
p <- 100
ci_method <- "quantile"
all_res <- list()
methods <- "zerosample2"
SNR <- 1
corr <- "exchangeable"
rho <- 0
alpha <- 0.2
niter <- 100

## Function
accross_lambda_res <- function(dat, ci_method, method) {

  truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))

  if (length(unique(abs(truth$truth))) <= 4) {
    group_cutoffs <- sort(unique(abs(truth$truth)))[-1]
    group_names <- sort(unique(abs(truth$truth)))
  } else if (mean(truth$truth == 0) > .33 ){
    group_cutoffs <- numeric(2)
    group_cutoffs[1] <- 0 + 1e-9
    group_cutoffs[2] <- quantile(abs(truth$truth[truth$truth != 0]), .5)
    print(group_cutoffs)
  } else {
    group_cutoffs <- quantile(abs(truth$truth), c(.333, .667))
  }

  group_names <- c("Small", "Moderate", "Large")

  lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
  lambda_min <- lambda_max * 0.001
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = 0.001)
  lambda_min <- cv_fit$lambda.min / lambda_max
  # if (!any(abs(lambda_min - lambda_seq) < 1e-4)) {
  #   lambda_seq <- sort(c(lambda_seq, lambda_min), decreasing = TRUE)
  # }

  res <- list()
  for (i in 1:length(lambda_seq)) {
    boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, method = method, max.iter = 1e8)
    res[[i]] <- ci.boot.ncvreg(boot_res, ci_method = ci_method) %>%
      dplyr::mutate(lambda = i) %>%
      dplyr::select(variable, lambda, lower, upper) %>%
      dplyr::left_join(truth) %>%
      dplyr::mutate(
        group = group_names[(1 + (abs(truth) >= group_cutoffs[1]) + (abs(truth) >= group_cutoffs[2]))]
      )
  }

  return(list("plot_data" = do.call("rbind", res), "lambda_min" = lambda_min))

}

for (i in 1:length(methods)) {

  method <- methods[i]
  avg_lambdas <- numeric(length(methods))

  ## Sparse 1
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    sparse_beta <- c(rep(-2, 1), rep(-1, 1), rep(-0.5, 3), rep(2, 1), rep(1, 1), rep(0.5, 3), rep(0, 90))
    # dat <- gen_data(n = n, p = p, beta = sparse_beta)
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = sparse_beta)
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[1]] <- do.call(rbind, tmp)
  avg_lambdas[1] <- mean(lambdas)


  ## Sparse 2
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    sparse_beta <- c(rnorm(30), rep(0, 70))
    # dat <- gen_data(n = n, p = p, beta = sparse_beta)
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = sparse_beta)
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[2]] <- do.call(rbind, tmp)
  avg_lambdas[2] <- mean(lambdas)

  ## Sparse 3
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    sparse_beta <- c(rnorm(50), rep(0, 50))
    # dat <- gen_data(n = n, p = p, beta = sparse_beta)
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = sparse_beta)
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[3]] <- do.call(rbind, tmp)
  avg_lambdas[3] <- mean(lambdas)

  ## Normal
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    # dat <- gen_data(n = n, p = p, beta = rnorm(p, mean = 0, sd = 1))
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rnorm(p, mean = 0, sd = 1))
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[4]] <- do.call(rbind, tmp)
  avg_lambdas[4] <- mean(lambdas)

  ## Laplace
  rt <- 2
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    # dat <- gen_data(n = n, p = p, beta = rlaplace(p, rate = rt))
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rlaplace(p, rate = rt))
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[5]] <- do.call(rbind, tmp)
  avg_lambdas[5] <- mean(lambdas)

  ## T
  set.seed(my_seed)
  tmp <- list()
  lambdas <- numeric(niter)
  for (j in 1:niter) {
    # dat <- gen_data(n = n, p = p, beta = rt(p, 3))
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rt(p, 3))
    res <- accross_lambda_res(dat, ci_method, method)
    tmp[[j]] <- res$plot_data %>% mutate(iter = j)
    lambdas[j] <- res$lambda_min
  }
  all_res[[6]] <- do.call(rbind, tmp)
  avg_lambdas[6] <- mean(lambdas)

  save(all_res, avg_lambdas, file = glue("./rds/distribution_of_beta_SNR{SNR}_{corr}_rho{rho*100}_alpha{alpha*100}_p{100}.rds"))
}


# n = 101 lm bootstrap
# - general problem
