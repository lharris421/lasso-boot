source("./scripts/setup/setup.R")

ns <- c(20, 30, 60)
p <- 30
rt <- 2
nboot <- 1000

per_var <- per_dataset <- list()
# method <- "bucketfill"
# methods <- c("selective_inference")
# methods <- "blp"
methods <- c("traditional", "zerosample1", "sample", "debiased", "acceptreject")
n_methods <- length(methods)

for (i in 1:length(ns)) {

  print(i)
  info <- list()
  info2 <- list()

  n <- ns[i]
  current_seed <- my_seed + n
  true_lambda <- (1 / n) * rt

  for (j in 1:100) {
    print(j)
    current_seed <- current_seed + j
    set.seed(current_seed)
    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data(n = n, p = p, beta = laplace_beta)
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    if (!(any(methods %in% c("selective_inference", "blp")))) {
      cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
      lambda <- cv_fit$lambda.min
      sigma2 <- cv_fit$cve[cv_fit$min]
    }

    coverage_df <- data.frame(laplace_beta, j)

    per_var_info <- list()
    per_dataset_info <- list()
    for (k in 1:n_methods) {

      set.seed(current_seed)

      start <- Sys.time()
      if (methods[k] == "selective_inference") {
        ### Selective Inference
        res <- selective_inference(dat, estimate_sigma = TRUE)
        if (!is.null(res)) {
          lam <- res$lambda
          ci <- dplyr::full_join(res$confidence_interval, data.frame(variable = names(dat$beta))) %>%
            mutate(method = methods[k])
        } else {
          lam <- NA
          ci <- data.frame(variable = names(dat$beta), lower = NA, upper = NA, estimate = NA)
        }
      } else if (methods[k] == "blp") {
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
        lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, quantiles = methods[k], max.iter = 1e6, lambda = lambda, sigma2 = sigma2)
        ci <- ci.boot.ncvreg(lassoboot, method = method)
        lam <- lambda
      }
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))


      per_var_info[[k]] <- inner_join(data.frame(truth = dat$beta, variable = names(dat$beta)), bind_cols(variable = ci$variable, ci$lower, ci$upper, ci$estimate, methods[k], j, n))
      per_dataset_info[[k]] <- bind_cols(elapsed, lam, methods[k], j, n)
    }

    info[[j]] <- do.call(rbind, per_var_info)
    info2[[j]] <- do.call(rbind, per_dataset_info)

  }

  per_var[[i]] <- do.call(rbind, info)
  per_dataset[[i]] <- do.call(rbind, info2)

}

per_var_all <- do.call(rbind, per_var)
colnames(per_var_all) <- c("truth", "variable", "lower", "upper", "estimate", "method", "group", "n")
per_dataset_all <- do.call(rbind, per_dataset)
colnames(per_dataset_all) <- c("time", "lambda", "method", "group", "n")


for (i in 1:length(methods)) {
  per_var <- per_var_all %>%
    filter(method == methods[i])
  per_dataset <- per_dataset_all %>%
    filter(method == methods[i])
  save(per_dataset, per_var, file = glue("./rds/laplace_{methods[i]}.rds"))
}

