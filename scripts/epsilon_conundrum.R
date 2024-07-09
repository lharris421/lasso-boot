source("./scripts/setup/setup.R")

method <- "lasso"
n_values <- 100
data_type <- "sparse"
SNR <- 1
alpha <- .2
p <- 100
modifier <- "debias"
enet_alpha <- 1
gamma <- NA
penalty <- method

args_list <- list(
  data = data_type, n = n_values, snr = SNR, lambda = "cv", modifier = modifier,
  method = method, nominal_coverage = (1-alpha) * 100, p = p,
  alpha = enet_alpha, gamma = gamma
)

rds_folder <- "/Users/loganharris/github/lasso-boot/rds"
check_parameters_existence(rds_folder, args_list, check_for = "existing", halt = TRUE)

example_it <- sample(1:100, 1)

## Data
res <- list()
current_seed <- my_seed
for (j in 1:100) {

  print(j)

  current_seed <- current_seed + j
  set.seed(current_seed)
  beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
  dat <- gen_data_snr(n = n_values, p = p, p1 = a, beta = beta, rho = 0, SNR = SNR)

  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = method)
  print(lambda <- cv_fit$lambda.min)
  print(sigma2 <- cv_fit$cve[cv_fit$min])

  set.seed(current_seed)

  tmp <- boot_ncvreg(
    dat$X, dat$y, lambda = lambda, sigma2 = sigma2, verbose = TRUE,
    nboot = nboot, max.iter = 1e6,
    penalty = method, alpha = enet_alpha, gamma = gamma, debias = TRUE)

  if (j == example_it) {
    example_res <- tmp
  }

  res[[j]] <- ci.boot_ncvreg(tmp, debias = TRUE) %>%
    mutate(submethod = method, method = penalty, group = j) %>%
    inner_join(data.frame(truth = dat$beta, variable = names(dat$beta)))

  res[[j]] %>%
    mutate(covered = truth >= lower & truth <= upper) %>%
    group_by(submethod) %>%
    summarise(coverage = mean(covered)) %>%
    print()
}

do.call(rbind, res) %>%
  mutate(covered = truth >= lower & truth <= upper) %>%
  group_by(submethod, abs(truth)) %>%
  summarise(coverage = mean(covered)) %>%
  print()

do.call(rbind, res) %>%
  mutate(covered = truth >= lower & truth <= upper) %>%
  group_by(submethod) %>%
  summarise(coverage = mean(covered)) %>%
  print()

confidence_interval <- res
example <- example_res
res_list <- list("confidence_interval" = confidence_interval, "example" = example)
save_objects(folder = rds_folder, res_list, args_list = args_list, overwrite = TRUE)

