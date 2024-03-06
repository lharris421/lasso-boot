source("./scripts/setup/setup.R")

## Parameters
data_type <- "abn"
alpha <- .2
modifier <- character()
modifier[1] <- NA
arg_list <- list(data = data_type,
                 n = 100,
                 p = 100,
                 snr = 1,
                 sd = ifelse(data_type == "normal", sd, NA),
                 rate = ifelse(data_type == "laplace", rt, NA),
                 a = ifelse(data_type == "abn", 1, NA),
                 b = ifelse(data_type == "abn", 1, NA),
                 correlation_structure = "exchangeable",
                 correlation = 0.99,
                 correlation_noise = ifelse(data_type == "abn", 0, NA),
                 method = c("ridge", "zerosample2"),
                 ci_method = "quantile",
                 nominal_coverage = alpha * 100,
                 modifier = modifier)

new_folder <- "/Users/loganharris/github/lasso-boot/new_rds"
check_parameters_existence(new_folder, arg_list, check_for = "existing", halt = TRUE)

## Selected example
selected_example <- sample(1:100, 1)

cis <- list()
examples <- list()
coverages <- list()
current_seed <- my_seed
for (i in 1:100) {

  print(i)
  current_seed <- current_seed + i
  set.seed(current_seed)
  dat <- gen_data_abn(n = arg_list$n, p = arg_list$p, a = arg_list$a, b = arg_list$b, rho = arg_list$correlation, SNR = arg_list$snr)

  for (j in 1:length(arg_list$method)) {
    set.seed(current_seed)
    if (i == 1) {
      cis[[j]] <- list()
    }
    if (arg_list$method[j] == "ridge") {
      ### Ridge
      ridge_cv <- cv_ridge(dat$X, dat$y)
      cis[[j]][[i]] <- cbind(confint(ridge_cv$fit, level = 1 - alpha, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min), method = "ridge") %>%
        data.frame() %>%
        mutate(group = i)
    } else {
      ### Lasso-boot sample
      lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, max.iter = 1e9, method = arg_list$method[j], nboot = nboot)
      cis[[j]][[i]] <- ci.boot.ncvreg(lassoboot, ci_method = arg_list$ci_method) %>%
        dplyr::mutate(method = arg_list$method[j], group = i)
    }

    coverages[[j]][[i]] <- cis[[j]][[i]] %>%
      mutate(truth = dat$beta, covered = lower <= truth & upper >= truth) %>%
      pull(covered) %>%
      mean()

    if (i == selected_example) {
      if (arg_list$method[j] == "ridge") {
        examples[[j]] <- cbind(confint(ridge_cv$fit, level = 1 - alpha, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))
      } else {
        examples[[j]] <- lassoboot
      }
    }

  }

}
names(cis) <- arg_list$method
names(examples) <- arg_list$method
names(coverages) <- arg_list$method

for (i in 1:length(arg_list$method)) {
  confidence_interval <- cis[[arg_list$method[i]]]
  example <- examples[[arg_list$method[i]]]
  alist <- list(data = data_type,
                   n = 100,
                   p = 100,
                   snr = 1,
                   sd = ifelse(data_type == "normal", sd, NA),
                   rate = ifelse(data_type == "laplace", rt, NA),
                   a = ifelse(data_type == "abn", 1, NA),
                   b = ifelse(data_type == "abn", 1, NA),
                   correlation_structure = "exchangeable",
                   correlation = 0.99,
                   correlation_noise = ifelse(data_type == "abn", 0, NA),
                   method = c("ridge", "zerosample2")[i],
                   ci_method = "quantile",
                   nominal_coverage = alpha * 100,
                   modifier = NA)
  save_objects(folder = new_folder, args_list = alist, confidence_interval, example, coverages, overwrite = TRUE)
}

