source("./scripts/setup/setup.R")

nboot <- 1000
simulations <- 100

alpha <- .2
args_list <- list(data = "laplace",
                    snr = 1,
                    n = 100,
                    p = 100,
                    method = "MCP", gamma = 3,
                    lambda = "across",
                    nominal_coverage = alpha * 100,
                    alpha = 1)

res <- list()
lambdas <- list()
true_lambdas <- list()
cov_lambdas <- list()
for (j in 1:length(args_list$n)) {

  n <- args_list$n[j]
  current_seed <- floor((my_seed + n) * alpha)

  lambda_mins <- numeric(simulations)
  lambda_truths <- numeric(simulations)
  nres <- list()
  for (k in 1:simulations) {

    current_seed <- current_seed + k
    set.seed(current_seed)
    laplace_beta <- rlaplace(args_list$p, rate = 1)
    dat <- gen_data_snr(n = n, p = args_list$p, p1 = args_list$p, beta = laplace_beta, rho = 0, SNR = args_list$snr)
    # dat <- gen_data_snr(n = n, p = args_list$p, p1 = args_list$p, beta = laplace_beta, rho = 0.8, corr = "autoregressive", SNR = args_list$snr)

    true_rate <- laplace_beta[1] / dat$beta[1]
    print(true_rate)
    true_lambda <- true_rate / n

    cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = args_list$method, lambda.min = 0.001)
    which_lambdas <- ceiling(((1:100)[1:100 %% 10 == 0 | 1:100 == 1] / 100) * length(cv_fit$lambda))
    lambda_seq <- cv_fit$lambda[which_lambdas]
    lambda_max <- max(cv_fit$lambda)
    lambda_mins[k] <- cv_fit$lambda.min / lambda_max
    lambda_truths[k] <- true_lambda / lambda_max
    print(k)
    print(lambda_mins[k])

    # if (!(cv_fit$lambda.min %in% lambda_seq)) {
    #   lambda_seq <- c(lambda_seq, cv_fit$lambda.min)
    # }

    pre_lambda_res <- list()
    for (i in 1:length(lambda_seq)) {
      set.seed(current_seed)
      boot_res <- boot_ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, method = args_list$method, lambda.min = 0.001, max.iter = 1e8)
      pre_lambda_res[[i]] <- ci.boot_ncvreg(boot_res, alpha = alpha) %>%
        dplyr::mutate(lambda_ind = i, n = n, group = k, lambda = lambda_seq[i], lambda_max = lambda_max)
    }

    truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))
    plot_data <- do.call("rbind", pre_lambda_res) %>%
      left_join(truth)

    plot_data %>%
      mutate(covered = truth >= lower & truth <= upper) %>%
      group_by(lambda_ind, method) %>%
      summarise(coverage = mean(covered)) %>%
      print()

    nres[[k]] <- plot_data

  }

  lambdas[[j]] <- lambda_mins
  true_lambdas[[j]] <- lambda_truths
  res[[j]] <- do.call(rbind, nres)

}

# args_list$modifier <- "correlation"
res_list <- list("res" = res, "lambdas" = lambdas, "true_lambdas" = true_lambdas)
save_objects(folder = rds_path, res_list, args_list = args_list, overwrite = TRUE, save_method = "rds")
