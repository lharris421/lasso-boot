source("./scripts/setup/setup.R")

## Data generating scenarios
# data_type <- "abn"
# a <- 5
# b <- 2
# rho <- .5
# rho.noise <- rho - .3
# corr <- "exchangeable"

# data_type <- "abn"
# a <- 5
# b <- 2
# rho <- .8
# rho.noise <- rho - .3
# corr <- "autoregressive"

# data_type <- "laplace"
# rt <- 2
# #rt <- 10
# corr <- "exchangeable"
# rho <- 0

# data_type <- "laplace"
# rt <- 2
# corr <- "autoregressive"
# rho <- .7

# data_type <- "normal"
# sd <- 1
# corr <- "exchangeable"
# rho <- 0

p <- 60
ns <- p * c(0.75, 1, 4)
nboot <- 1000
simulations <- 100
alpha <- .2
SNR <- 1

# method <- c("zerosample1", "zerosample2", "sample", "debiased", "acceptreject", "traditional", "blp", "selective_inference")
# method <- c("zerosample1", "zerosample2", "sample", "debiased", "acceptreject", "traditional")
# method <- c("sample", "debiased", "traditional")
# method <- c("fullconditional")
method <- c("zerosample2")
# method <- c("truncatedzs2")
# method <- c("selective_inference")
nboot <- ifelse(all(method == "fullconditional"), 1, nboot)
ci_method <- ifelse(all(method == "fullconditional"), "identity", ci_method)
n_method <- length(method)

per_var <- per_dataset <- list()
for (i in 1:length(ns)) {

  print(i)
  info <- list()
  info2 <- list()

  n <- ns[i]
  current_seed <- floor((my_seed + n) * alpha)

  for (j in 1:simulations) {
    print(j)
    current_seed <- current_seed + j
    set.seed(current_seed)

    if (data_type == "laplace") {
      true_lambda <- (1 / n) * rt
      laplace_beta <- rlaplace(p, rate = rt)
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "abn") {
      dat <- gen_data_abn(n = n, p = p, a = a, b = b, rho = rho, rho.noise = rho.noise, noise = corr, SNR = SNR)
    } else if (data_type == "normal") {
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rnorm(p, sd = sd), corr = corr, rho = rho, SNR = SNR)
    }
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    ## Only needed for running true lambda / sigma
    # lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
    # lambda_min <- (true_lambda / lambda_max) * .999
    # lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))

    ## Estimate lambda / sigma2
    if (any(!(method %in% c("selective_inference", "blp")))) {
      cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", max.iter = 1e8)
      # cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = lambda_min, max.iter = 1e8)
      lambda <- cv_fit$lambda.min
      sigma2 <- cv_fit$cve[cv_fit$min]
    }
    # print(sigma2)

    ## True lambda / CVE estimate sigma2
    # lambda <- true_lambda
    # ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
    # l <- floor(ind)
    # r <- ceiling(ind)
    # w <- ind %% 1
    # sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]

    ## True lambda / sigma2
    # lambda <- true_lambda
    # sigma2 <- drop(crossprod(dat$beta)) / SNR
    # print(sigma2)

    coverage_df <- data.frame(dat$beta, j)

    per_var_info <- list()
    per_dataset_info <- list()
    for (k in 1:n_method) {

      set.seed(current_seed)

      start <- Sys.time()
      if (method[k] == "selective_inference") {
        ### Selective Inference
        res <- selective_inference(dat, estimate_sigma = TRUE)
        if (!is.null(res)) {
          lam <- res$lambda
          ci <- dplyr::full_join(res$confidence_interval, data.frame(variable = names(dat$beta))) %>%
            mutate(method = method[k])
        } else {
          lam <- NA
          ci <- data.frame(variable = names(dat$beta), lower = NA, upper = NA, estimate = NA)
        }
      } else if (method[k] == "blp") {
        ### HDI - Across a range of lambda values
        res <- blp(dat)
        if (!is.null(res)) {
          lam <- res$lambda
          ci <- res$confidence_interval
        } else {
          lam <- NA
          ci <- data.frame(variable = names(dat$beta), lower = NA, upper = NA, estimate = NA)
        }
      } else {
        # lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, method = method[k], max.iter = 1e6, lambda = lambda, sigma2 = sigma2, lambda.min = lambda_min)
        lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, method = method[k], max.iter = 1e8, lambda = lambda, sigma2 = sigma2) ## add alpha if running full conditional
        ci <- ci.boot.ncvreg(lassoboot, ci_method = ci_method, alpha = alpha, original_data = dat)
        lam <- lambda

      }
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))


      per_var_info[[k]] <- inner_join(data.frame(truth = dat$beta, variable = names(dat$beta)), bind_cols(variable = ci$variable, lower = ci$lower, upper = ci$upper, ci$estimate, method[k], j, n))
      per_dataset_info[[k]] <- bind_cols(elapsed, lam, method[k], j, n)
      print(per_var_info[[k]] %>% mutate(covered = lower <= truth & upper >= truth) %>% pull(covered) %>% mean(na.rm = TRUE))
    }

    info[[j]] <- do.call(rbind, per_var_info)
    info2[[j]] <- do.call(rbind, per_dataset_info)

  }

  per_var[[i]] <- do.call(rbind, info)
  per_dataset[[i]] <- do.call(rbind, info2)
  print(per_var[[i]] %>% mutate(covered = lower <= truth & upper >= truth) %>% pull(covered) %>% mean(na.rm = TRUE))

}

per_var_all <- do.call(rbind, per_var)
colnames(per_var_all) <- c("truth", "variable", "lower", "upper", "estimate", "method", "group", "n")
per_dataset_all <- do.call(rbind, per_dataset)
colnames(per_dataset_all) <- c("time", "lambda", "method", "group", "n")


for (i in 1:length(method)) {
  per_var <- per_var_all %>%
    filter(method == method[i])
  per_dataset <- per_dataset_all %>%
    filter(method == method[i])
  ad_inf <- switch(
    data_type,
    laplace = rt,
    abn = a,
    normal = sd
  )
  save(per_dataset, per_var, file = glue("./rds/{data_type}({ad_inf})_SNR{SNR}_{corr}_rho{rho*100}_{method[i]}_alpha{alpha*100}_p{p}.rds"))
}
