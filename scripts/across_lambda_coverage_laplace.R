source("./scripts/setup/setup.R")

nboot <- 1000
simulations <- 100

rt <- 2
# ns <- c(30, 60, 120)
p <- 100
ns <- p * c(1)
method <- "zerosample2"
ci_method <- "quantile"
alpha <- .2
SNR <- 1
corr <- "exchangeable"
rho <- 0

res <- list()
lambdas <- list()
for (j in 1:length(ns)) {

  n <- ns[j]
  true_lambda <- (1 / n) * rt
  current_seed <- floor((my_seed + n) * alpha)

  lambda_mins <- numeric(simulations)
  nres <- list()
  for (k in 1:simulations) {

    current_seed <- current_seed + k
    set.seed(current_seed)
    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)

    lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
    lambda_min <- lambda_max * 0.001
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

    cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = 0.001)
    lambda_mins[k] <- cv_fit$lambda.min / lambda_max
    print(k)
    print(lambda_mins[k])

    # if (!any(abs(lambda_min - lambda_seq) < 1e-4)) {
    #   lambda_seq <- sort(c(lambda_seq, lambda_min), decreasing = TRUE)
    # }

    pre_lambda_res <- list()
    for (i in 1:length(lambda_seq)) {
      set.seed(current_seed)
      boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, method = method, lambda.min = 0.001, max.iter = 1e8)
      pre_lambda_res[[i]] <- ci.boot.ncvreg(boot_res, ci_method = ci_method) %>%
        dplyr::mutate(lambda_ind = i, n = n, group = k)
    }

    truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))
    plot_data <- do.call("rbind", pre_lambda_res) %>%
      left_join(truth)

    nres[[k]] <- plot_data

  }

  lambdas[[j]] <- lambda_mins
  res[[j]] <- do.call(rbind, nres)

}

save(res, lambdas, file = glue("./rds/across_lambda_laplace({rt})_SNR{SNR}_{corr}_rho{rho*100}_{method}_alpha{alpha*100}_p{100}.rds"))
