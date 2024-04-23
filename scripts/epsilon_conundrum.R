source("./scripts/setup/setup.R")

methods <- c("zerosample2", "traditional")
n_values <- 100
data_type <- "sparse"
SNR <- 1
alpha <- .2
p <- 100
ci_method <- "quantile"

args_list <- list(data = data_type, n = n_values, snr = SNR, lambda = "cv",
                                method = methods,
                                ci_method = ci_method, nominal_coverage = alpha * 100, p = p)

rds_folder <- "/Users/loganharris/github/lasso-boot/rds"
check_parameters_existence(rds_folder, args_list, check_for = "existing", halt = TRUE)

example_it <- sample(1:100, 1)

## Data
res <- list()
example_res <- list()
current_seed <- my_seed
for (j in 1:100) {
  print(j)

  current_seed <- current_seed + j
  set.seed(current_seed)
  beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
  dat <- gen_data_snr(n = n_values, p = p, p1 = a, beta = beta, rho = 0, SNR = SNR)

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  print(lambda <- cv_fit$lambda.min)
  print(sigma2 <- cv_fit$cve[cv_fit$min])

  for (i in 1:length(methods)) {
    set.seed(current_seed)
    if (j == 1) {res[[i]] <- list()}
    tmp <- boot.ncvreg(dat$X, dat$y, lambda = lambda, sigma2 = sigma2, verbose = TRUE, nboot = nboot, method = methods[i], max.iter = 1e6)

    if (j == example_it) {
      example_res[[i]] <- tmp
    }

    res[[i]][[j]] <- ci.boot.ncvreg(tmp, ci_method = ci_method, original_data = dat) %>%
      mutate(method = methods[i], truth = dat$beta, group = j)

    res[[i]][[j]] %>%
      mutate(covered = truth >= lower & truth <= upper) %>%
      pull(covered) %>%
      mean() %>%
      print()
  }
}
names(res) <- methods
names(example_res) <- methods

do.call(rbind, res[[i]]) %>%
  mutate(covered = truth >= lower & truth <= upper) %>%
  pull(covered) %>%
  mean() %>%
  print()


for (i in 1:length(methods)) {
    confidence_interval <- res[[methods[i]]]
    example <- example_res[[methods[i]]]
    res_list <- list("confidence_interval" = confidence_interval, "example" = example)
    args_list <- list(data = data_type, n = n_values, snr = SNR, lambda = "cv",
                      method = methods[i],
                      ci_method = ci_method, nominal_coverage = alpha * 100, p = p)
    save_objects(folder = rds_folder, res_list, args_list = args_list, overwrite = TRUE)
}

