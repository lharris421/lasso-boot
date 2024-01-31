source("./scripts/setup/setup.R")

ns <- c(20, 30, 60)
p <- 30
rt <- 2
nboot <- 1000

res_width <- res_coverage <- res_time <- res_lambda <- list()
all_info <- list()
# method <- "bucketfill"

for (i in 1:length(ns)) {

  print(i)
  lambdas <- matrix(nrow = 100, ncol = n_methods)
  times <- matrix(nrow = 100, ncol = n_methods)
  widths <- matrix(nrow = 100, ncol = n_methods)
  coverages <- list()
  info <- list()

  n <- ns[i]
  true_lambda <- (1 / n) * rt

  for (j in 1:100) {
    print(j)
    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data(n = n, p = p, beta = laplace_beta)
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    coverage_df <- data.frame(laplace_beta, j)

    methods_info <- list()
    for (k in 1:n_methods) {
      quantiles <- methods[k]
      start <- Sys.time()
      lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, quantiles = quantiles, max.iter = 1e6)
      times[j,k] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      lassoboot_ci <- ci.boot.ncvreg(lassoboot, method = method)
      lassoboot_coverage <- lassoboot_ci$lower <= laplace_beta & laplace_beta <= lassoboot_ci$upper
      lambdas[j,k] <- lassoboot$lambda
      widths[j,k] <- median(lassoboot_ci$upper - lassoboot_ci$lower)
      coverage_df <- bind_cols(coverage_df, lassoboot_coverage)
      tmp <- bind_cols(laplace_beta, lassoboot_ci$lower, lassoboot_ci$upper, lassoboot$estimates, quantiles, j)
      colnames(tmp) <- c("truth", "lower", "upper", "estimate", "method", "group")
      methods_info[[k]] <- tmp
    }

    colnames(coverage_df) <- c("truth", "group", methods)
    coverages[[j]] <- coverage_df
    info[[j]] <- do.call(rbind, methods_info)

  }

  res_time[[i]] <- times
  res_coverage[[i]] <- coverages
  res_lambda[[i]] <- lambdas
  res_width[[i]] <- widths
  all_info[[i]] <- do.call(rbind, info)

}

save(res_time, res_coverage, res_lambda, res_width, ns, all_info, file = glue("./rds/lassoboot_comparison_laplace_100_all_{method}_1000.rds"))
