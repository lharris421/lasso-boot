source("./scripts/setup/setup.R")
library(tictoc)
## Data arguments
data_type <- "laplace"
rt <- 2
corr <- "exchangeable"
rho <- 0
# rho.noise <- 0
# a <- 5
# b <- 2
# sd <- 1
p <- 100
ns <- p * c(0.5, 1, 4)
nboot <- 1000
simulations <- 100
alpha <- .2
SNR <- 1
modifier <- NA

methods <- c("zerosample2")
n_methods <- length(methods)
ci_method <- "quantile"

args_list <- list(data = data_type,
                  n = ns,
                  snr = SNR,
                  sd = ifelse(data_type == "normal", sd, NA),
                  rate = ifelse(data_type == "laplace", rt, NA),
                  a = ifelse(data_type == "abn", a, NA),
                  b = ifelse(data_type == "abn", b, NA),
                  correlation_structure = corr,
                  correlation = rho * 100,
                  correlation_noise = ifelse(data_type == "abn", rho.noise * 100, NA),
                  method = methods,
                  ci_method = ci_method,
                  nominal_coverage = alpha * 100,
                  lambda = "cv_1se",
                  # modifier = modifier,
                  p = p)

# check_parameters_existence(rds_folder, args_list, check_for = "existing", halt = TRUE)

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
    } else if (data_type == "orthogonal") {
      true_lambda <- (1 / n) * rt
      laplace_beta <- rlaplace(p, rate = rt)
      dat <- genOrthoSNR(n = n, p = p, beta = laplace_beta, SNR = SNR)
    }
    # dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    ## Only needed for running true lambda / sigma
    if (!is.na(modifier)) {
      lambda_max <- max(ncvreg:::find_thresh(std(dat$X), dat$y))
      lambda_min <- (true_lambda / lambda_max) * .999
    }

    ## Estimate lambda / sigma2
    if (any(!(methods %in% c("selectiveinference", "blp")))) {
      cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", lambda.min = ifelse(is.na(modifier), ifelse(n > p, 0.001, 0.05), lambda_min), max.iter = 1e8)

      ## True lambda / CVE estimate sigma2
      if (!is.na(modifier)) {
        lambda <- true_lambda
        if (modifier == "tl") {
          ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
          l <- floor(ind)
          r <- ceiling(ind)
          w <- ind %% 1
          sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
        } else if (modifier == "tls") {
          sigma2 <- drop(crossprod(dat$beta)) / SNR
        }
      } else {
        lambda <- cv_fit$lambda.min
        # lambda <- cv_fit$lambda[min(which(cv_fit$lambda <= (cv_fit$lambda.min + cv_fit$cvse[cv_fit$min])))]
        sigma2 <- cv_fit$cve[cv_fit$min]
      }

    }

    coverage_df <- data.frame(dat$beta, j)

    per_var_info <- list()
    per_dataset_info <- list()
    for (k in 1:n_methods) {

      set.seed(current_seed)
      print(methods[k])
      start <- Sys.time()
      if (methods[k] == "selectiveinference") {
        ### Selective Inference
        dat$X <- ncvreg::std(dat$X)
        res <- selective_inference(dat, estimate_sigma = FALSE)
        if (!is.null(res)) {
          lam <- res$lambda
          ci <- dplyr::full_join(res$confidence_interval, data.frame(variable = names(dat$beta))) %>%
            mutate(method = methods[k])
        } else {
          lam <- NA
          ci <- data.frame(variable = names(dat$beta), lower = NA, upper = NA, estimate = NA)
        }
      } else if (methods[k] == "blp") {
        ### HDI
        res <- blp(dat)
        if (!is.null(res)) {
          lam <- res$lambda
          ci <- res$confidence_interval
        } else {
          lam <- NA
          ci <- data.frame(variable = names(dat$beta), lower = NA, upper = NA, estimate = NA)
        }
      } else {

        if (is.na(modifier)) {
          lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = TRUE, nboot = nboot, method = methods[k], max.iter = 1e8, lambda = lambda, sigma2 = sigma2) ## add alpha if running full conditional
        } else {
          lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, method = methods[k], max.iter = 1e8, lambda = lambda, sigma2 = sigma2, lambda.min = lambda_min)
        }
        ci <- ci.boot.ncvreg(lassoboot, ci_method = ci_method, alpha = alpha, original_data = dat)
        lam <- lambda

      }
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

      lam <- ifelse(is.null(lam), NA, lam)
      print(methods[k])
      per_var_info[[k]] <- inner_join(data.frame(truth = dat$beta, variable = names(dat$beta)), bind_cols(variable = ci$variable, lower = ci$lower, upper = ci$upper, ci$estimate, methods[k], j, n))
      per_dataset_info[[k]] <- bind_cols(elapsed, lam, methods[k], j, n)
      print(per_dataset_info[[k]])
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


per_var_all %>%
  mutate(covered = truth >= lower & truth <= upper) %>%
  group_by(n) %>%
  summarise(coverage = mean(covered))

for (i in 1:length(methods)) {
  for (j in 1:length(ns)) {
    per_var_n <- per_var_all %>%
      filter(method == methods[i] & n == ns[j])
    per_dataset_n <- per_dataset_all %>%
      filter(method == methods[i] & n == ns[j])
    args_list <- list(data = data_type,
         n = ns[j],
         snr = SNR,
         sd = ifelse(data_type == "normal", sd, NA),
         rate = ifelse(data_type == "laplace", rt, NA),
         a = ifelse(data_type == "abn", a, NA),
         b = ifelse(data_type == "abn", b, NA),
         # correlation_structure = corr,
         # correlation = rho * 100,
         correlation_noise = ifelse(data_type == "abn", rho.noise * 100, NA),
         method = methods[i],
         ci_method = ci_method,
         nominal_coverage = alpha * 100,
         lambda = "cv",
         # modifier = modifier,
         p = p)
    save_objects(folder = rds_folder, args_list = args_list, overwrite = TRUE, per_var_n, per_dataset_n)
  }
}

