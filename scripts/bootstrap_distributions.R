source("./scripts/setup/setup.R")
library(tictoc)

## Data arguments
data_type <- "laplace"
rt <- 2
p <- 100
ns <- 100
nboot <- 1000
alpha <- .2
SNR <- 1

methods <- c("zerosample2")
n_methods <- length(methods)
ci_method <- "quantile"

args_list <- list(data = data_type,
                  n = ns,
                  snr = SNR,
                  rate = rt,
                  method = methods,
                  ci_method = ci_method,
                  nominal_coverage = alpha * 100,
                  lambda = "cv",
                  p = p)

n <- ns
current_seed <- floor((my_seed + n) * alpha)
true_lambda <- (1 / n) * rt
laplace_beta <- rlaplace(p, rate = rt)
dat <- gen_data_snr(n = n, p = p, p1 = p, beta = laplace_beta, SNR = SNR)
truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)
res <- ncvreg::boot.ncvreg(dat$X, dat$y, nboot = 1000)

save_objects(folder = rds_path, res, truth_df,  args_list = args_list, overwrite = FALSE, save_method = "rda")
