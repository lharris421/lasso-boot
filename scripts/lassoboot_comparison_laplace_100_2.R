source("./scripts/setup/setup.R")

ns <- c(20, 30, 60)
p <- 30
rt <- 2
nboot <- 10000

method <- "quantile"

# res_width <- res_coverage <- res_time <- res_lambda <- matrix(nrow = 3, ncol = 3)
res_width <- res_coverage <- res_time <- res_lambda <- list()
methods <- c("mode", "sample", "zs", "disturbed")
n_methods <- 4

for (i in 1:length(ns)) {

  print(i)
  lambdas <- matrix(nrow = 100, ncol = n_methods)
  times <- matrix(nrow = 100, ncol = n_methods)
  widths <- matrix(nrow = 100, ncol = n_methods)
  coverages <- list()

  n <- ns[i]
  true_lambda <- (1 / n) * rt

  for (j in 1:100) {
    print(j)
    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data(n = n, p = p, beta = laplace_beta)
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    coverage_df <- data.frame(laplace_beta, j)

    for (k in 1:n_methods) {
      quantiles <- methods[k]
      start <- Sys.time()
      lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, quantiles = quantiles)
      times[j,k] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      lassoboot_ci <- ci.boot.ncvreg(lassoboot, method = method, original_data = dat)
      lassoboot_coverage <- lassoboot_ci$lower <= laplace_beta & laplace_beta <= lassoboot_ci$upper
      lambdas[j,k] <- lassoboot$lambda
      widths[j,k] <- median(lassoboot_ci$upper - lassoboot_ci$lower)
      coverage_df <- bind_cols(coverage_df, lassoboot_coverage)
    }

    colnames(coverage_df) <- c("truth", "group", methods)
    coverages[[j]] <- coverage_df

  }

  res_time[[i]] <- times
  res_coverage[[i]] <- coverages
  res_lambda[[i]] <- lambdas
  res_width[[i]] <- widths

}

save(res_time, res_coverage, res_lambda, res_width, ns, file = glue("./rds/lassoboot_comparison_laplace_100_all_{method}_2.rds"))
