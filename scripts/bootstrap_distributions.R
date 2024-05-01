source("./scripts/setup/setup.R")
library(tictoc)

## Data arguments
data_type <- "laplace"
p <- 100
ns <- 100
nboot <- 1000
alpha <- .2
SNR <- 1

methods <- c("sample")
n_methods <- length(methods)
ci_method <- "quantile"

args_list <- list(data = data_type,
                  n = ns,
                  snr = SNR,
                  method = methods,
                  ci_method = ci_method,
                  nominal_coverage = alpha * 100,
                  lambda = "cv",
                  p = p)

n <- ns
current_seed <- floor((my_seed + n) * alpha) + 1
set.seed(current_seed)
laplace_beta <- rlaplace(p, rate = 1)
dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, SNR = SNR)
truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)
res <- ncvreg::boot.ncvreg(dat$X, dat$y, nboot = 1000, method = args_list$method)
res_list <- list("res" = res, "truth_df" = truth_df)

args_list$data <- "laplace-single"
save_objects(folder = rds_path, res_list, args_list = args_list, overwrite = TRUE, save_method = "rds")
