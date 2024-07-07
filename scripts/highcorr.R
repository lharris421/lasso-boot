source("./scripts/setup/setup.R")

## Parameters
data_type <- "abn"
penalty <- method <- "ridge"
alpha <- .2
modifier <- NA
enet_alpha <- NA
gamma <- NA
arg_list <- list(data = data_type,
                 n = 100,
                 p = 100,
                 snr = 1,
                 a = ifelse(data_type == "abn", 1, NA),
                 b = ifelse(data_type == "abn", 1, NA),
                 correlation_structure = "exchangeable",
                 correlation = 0.99,
                 correlation_noise = ifelse(data_type == "abn", 0, NA),
                 method = method,
                 nominal_coverage = (1-alpha) * 100,
                 alpha = enet_alpha,
                 gamma = gamma,
                 lambda = "cv",
                 modifier = modifier)

rds_folder <- "/Users/loganharris/github/lasso-boot/rds"
check_parameters_existence(rds_folder, arg_list, check_for = "existing", halt = TRUE)

## Selected example
selected_example <- sample(1:100, 1)

cis <- list()
coverages <- list()
current_seed <- my_seed
for (i in 1:100) {

  print(i)
  current_seed <- current_seed + i
  set.seed(current_seed)
  dat <- hdrm::gen_data_abn(n = arg_list$n, p = arg_list$p, a = arg_list$a, b = arg_list$b, rho = arg_list$correlation, SNR = arg_list$snr)

  set.seed(current_seed)

  if (method == "ridge") {
    ### Ridge
    ridge_cv <- cv_ridge(dat$X, dat$y)
    cis[[i]] <- cbind(confint(ridge_cv$fit, level = 1 - alpha, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min), method = "ridge")[-1,] %>%
      data.frame() %>%
      mutate(group = i)
    cis[[i]]$variable <- rownames(cis[[i]])
    colnames(cis[[i]]) <- c("lower", "upper", "estimate", "method", "group", "variable")
  } else {
    ### Lasso-boot sample
    lassoboot <- boot_ncvreg(
      dat$X, dat$y, verbose = FALSE, max.iter = 1e9, nboot = nboot,
      penalty = method, alpha = enet_alpha, gamma = gamma
    )
    cis[[i]] <- ci.boot_ncvreg(lassoboot, alpha = alpha) %>%
      dplyr::mutate(submethod = method, method = penalty, group = i)
  }

  coverages[[i]] <- cis[[i]] %>%
    left_join(data.frame(truth = dat$beta, variable = names(dat$beta))) %>%
    mutate(covered = lower <= truth & upper >= truth, submethod = method) %>%
    group_by(submethod) %>%
    summarise(coverage = mean(covered))
  print(coverages[[i]])

  if (i == selected_example) {
    if (method == "ridge") {
      example <- cbind(confint(ridge_cv$fit, level = 1 - alpha, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))
      colnames(example) <- c("lower", "upper", "estimate")
    } else {
      example <- lassoboot
    }
  }

}

res_list <- list(confidence_interval = cis, example = example, coverages = coverages)
save_objects(folder = rds_folder, args_list = arg_list, res_list, overwrite = TRUE)

