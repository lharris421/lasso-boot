source("./scripts/setup/setup.R")

nboot <- 1000
simulations <- 100

alpha <- .2
args_list <- list(data = "laplace",
                    rate = 2,
                    snr = 1,
                    n = 100,
                    p = 100,
                    correlation_structure = "exchangeable",
                    correlation = 0,
                    correlation_noise = NA,
                    method = "debiased",
                    ci_method = "mvn",
                    lambda = "across",
                    nominal_coverage = alpha * 100)

res <- list()
lambdas <- list()
for (j in 1:length(args_list$n)) {

  n <- args_list$n[j]
  true_lambda <- (1 / n) * args_list$rate
  current_seed <- floor((my_seed + n) * alpha)

  lambda_mins <- numeric(simulations)
  nres <- list()
  for (k in 1:simulations) {

    current_seed <- current_seed + k
    set.seed(current_seed)
    laplace_beta <- rlaplace(args_list$p, rate = args_list$rate)
    dat <- gen_data_snr(n = n, p = args_list$p, p1 = args_list$p, beta = laplace_beta, corr = args_list$correlation_structure, rho = args_list$correlation, SNR = args_list$snr)

    lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
    lambda_min <- lambda_max * 0.001
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

    cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = 0.001)
    lambda_mins[k] <- cv_fit$lambda.min / lambda_max
    print(k)
    print(lambda_mins[k])

    lambda_seq <- c(lambda_seq, cv_fit$lambda.min)

    pre_lambda_res <- list()
    for (i in 1:length(lambda_seq)) {
      set.seed(current_seed)
      boot_res <- boot.ncvreg(X = dat$X, y = dat$y, lambda = lambda_seq[i], nboot = nboot, method = args_list$method, lambda.min = 0.001, max.iter = 1e8)
      pre_lambda_res[[i]] <- ci.boot.ncvreg(boot_res, ci_method = args_list$ci_method) %>%
        dplyr::mutate(lambda_ind = i, n = n, group = k)
    }

    truth <- data.frame(variable = names(dat$beta), truth = as.numeric(dat$beta))
    plot_data <- do.call("rbind", pre_lambda_res) %>%
      left_join(truth)

    plot_data %>%
      mutate(covered = truth >= lower & truth <= upper) %>%
      group_by(lambda_ind) %>%
      summarise(coverage = mean(covered)) %>%
      print()

    nres[[k]] <- plot_data

  }

  lambdas[[j]] <- lambda_mins
  res[[j]] <- do.call(rbind, nres)

}

save_objects(folder = rds_folder, args_list = args_list, overwrite = FALSE, res, lambdas)
