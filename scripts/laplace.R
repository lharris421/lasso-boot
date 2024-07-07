source("./scripts/setup/setup.R")
library(tictoc)
## Data arguments
data_type <- "laplace"
corr <- "autoregressive"
rho <- 0.8
p <- 100
# ns <- p * c(0.5, 1, 4)
ns <- p * 1
nboot <- 1000
simulations <- 100
alpha <- 0.1
SNR <- 1
modifier <- NA
enet_alpha <- 1
gamma <- NA

penalty <- method <- "lasso"
n_methods <- length(methods)

args_list <- list(data = data_type,
                  n = ns,
                  snr = SNR,
                  correlation_structure = ifelse(rho > 0, corr, NA),
                  correlation = ifelse(!is.null(corr), rho * 100, NA),
                  method = method,
                  nominal_coverage = (1-alpha) * 100,
                  lambda = "cv",
                  p = p,
                  alpha = enet_alpha,
                  gamma = gamma,
                  modifier = modifier)

# check_parameters_existence(rds_path, args_list, check_for = "existing", halt = FALSE)
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

      laplace_beta <- rlaplace(p, rate = 1)
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, corr = corr, rho = rho, SNR = SNR)

      true_rate <- sqrt(2*p)
      true_lambda <- true_rate / n

    } else if (data_type == "normal") {
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rnorm(p, sd = 2), corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "t") {
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rt(p, df = 4), corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "orthogonal") {
      dat <- hdrm::genOrtho(n = n, p = p)
    } else if (data_type == "uniform") {
      unif_beta <- runif(p, -1, 1)
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = unif_beta, corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "beta") {
      beta_beta <- rbeta(p, .1, .1) - .5
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = beta_beta, corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "sparse 1") {
      betas <- c(rep(c(rep(0.5, 3), 1, 2), 2) * c(rep(1, 5), rep(-1, 5)), rep(0, 90))
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = betas, corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "sparse 2") {
      betas <- c(rnorm(30), rep(0, 70))
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = betas, corr = corr, rho = rho, SNR = SNR)
    } else if (data_type == "sparse 3") {
      betas <- c(rnorm(50), rep(0, 50))
      dat <- gen_data_snr(n = n, p = p, p1 = p, beta = betas, corr = corr, rho = rho, SNR = SNR)
    }

    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)


    ## Estimate lambda / sigma2
    if (method == "lasso") {

      ## Only needed for running true lambda / sigma
      if (!is.na(modifier)) {
        lambda_max <- max(ncvreg(dat$X, dat$y, penalty = "lasso")$lambda)
        lambda_min <- (true_lambda / lambda_max) * .999
      }

      cv_fit <- cv.ncvreg(
        dat$X, dat$y, penalty = method, lambda.min = ifelse(is.na(modifier), ifelse(n > p, 0.001, 0.05), lambda_min),
        max.iter = 1e8, alpha = enet_alpha, gamma = gamma
      )

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
          sigma2 <- 1
        }
      } else {
        lambda <- cv_fit$lambda.min
        sigma2 <- cv_fit$cve[cv_fit$min]
      }

    }

    per_var_info <- list()
    per_dataset_info <- list()

    set.seed(current_seed)
    start <- Sys.time()
    if (method == "selectiveinference") {
      ### Selective Inference
      dat$X <- ncvreg::std(dat$X)
      res <- selective_inference(dat, estimate_sigma = FALSE)
      if (!is.null(res)) {
        lam <- res$lambda
        ci <- dplyr::full_join(res$confidence_interval, data.frame(variable = names(dat$beta))) %>%
          select(-method) %>% select(-scale)
      } else {
        lam <- NA
        ci <- data.frame(lower = NA, upper = NA, variable = names(dat$beta), estimate = NA)
      }
    } else if (method == "blp") {
      ### HDI
      res <- blp(dat)
      if (!is.null(res)) {
        print("Res is not null")
        lam <- res$lambda
        ci <- res$confidence_interval %>% select(-method)
      } else {
        print("Res is null")
        lam <- NA
        ci <- data.frame(lower = NA, upper = NA, variable = names(dat$beta), estimate = NA)
      }
    } else {

      if (is.na(modifier)) {
        lassoboot <- boot_ncvreg(
          dat$X, dat$y, penalty = method, verbose = TRUE, nboot = nboot,
          max.iter = 1e8, lambda = lambda, sigma2 = sigma2, debias = TRUE,
          alpha = enet_alpha, gamma = gamma)
      } else {
        lassoboot <- boot_ncvreg(
          dat$X, dat$y, penalty = method, verbose = FALSE, nboot = nboot, debias = TRUE,
          max.iter = 1e8, lambda = lambda, sigma2 = sigma2, lambda.min = lambda_min,
          alpha = enet_alpha, gamma = gamma)
      }
      ci <- ci.boot_ncvreg(lassoboot, alpha = alpha, debias = TRUE)
      lam <- lambda

    }
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    lam <- ifelse(is.null(lam), NA, lam)
    if (method %in% c("blp", "selectiveinference")) {
      per_var_info <- left_join(data.frame(truth = dat$beta, variable = names(dat$beta)), bind_cols(ci, method, j, n))
      colnames(per_var_info)[6:8] <- c("method", "group", "n")
    } else {
      per_var_info <- left_join(data.frame(truth = dat$beta, variable = names(dat$beta)), bind_cols(ci, j, n))
    }

    per_dataset_info <- bind_cols(elapsed, lam, method, j, n)
    print(per_var_info %>%
            mutate(covered = lower <= truth & upper >= truth) %>%
            group_by(method) %>%
            summarise(coverage = mean(covered, na.rm = TRUE))
    )

    info[[j]] <- per_var_info
    info2[[j]] <- per_dataset_info

  }

  per_var[[i]] <- do.call(rbind, info)
  per_dataset[[i]] <- do.call(rbind, info2)
  print(per_var[[i]] %>%
          mutate(covered = lower <= truth & upper >= truth) %>%
          group_by(method) %>%
          summarise(coverage = mean(covered, na.rm = TRUE))
  )

}

per_var_all <- do.call(rbind, per_var)
colnames(per_var_all) <- c("truth", "variable", "estimate", "lower", "upper", "method", "group", "n")
# colnames(per_var_all) <- c("truth", "variable", "lower", "upper", "estimate", "method", "group", "n")
per_dataset_all <- do.call(rbind, per_dataset)
colnames(per_dataset_all) <- c("time", "lambda", "method", "group", "n")

per_var_all <- per_var_all %>%
  select(variable, truth, lower, upper, estimate, method, group, n)


per_var_all %>%
  mutate(covered = truth >= lower & truth <= upper) %>%
  group_by(n, method) %>%
  summarise(coverage = mean(covered, na.rm = TRUE))

for (j in 1:length(ns)) {
  per_var_n <- per_var_all %>%
    rename(submethod = method) %>%
    mutate(method = "lasso") %>%
    # mutate(submethod = method, method = penalty) %>%
    filter(n == ns[j])
  per_dataset_n <- per_dataset_all %>%
    filter(n == ns[j])

  tmp_args <- args_list
  tmp_args$n <- ns[j]
  tmp_args$modifier <- "debias_z"

  res_list <- list("per_var_n" = per_var_n, "per_dataset_n" = per_dataset_n)
  save_objects(folder = rds_path, res_list, args_list = tmp_args, overwrite = TRUE, save_method = "rds")
}

