source("./scripts/setup/setup.R")

n <- 50
p <- 50
quantiles <- "zerosample"
method <- "quantile"
dist <- "laplace"

coverages <- list()
for (j in 1:100) {

  if (dist == "laplace") {
    beta <- rlaplace(p, rate = 2)
  } else if (dist == "normal") {
    beta <- rnorm(p, 0, 1)
  }

  dat <- gen_data(n = n, p = p, beta = beta)
  truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

  ### Lasso-boot
  lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = 1000, quantiles = quantiles)
  lassoboot_ci <- ci.boot.ncvreg(lassoboot, alpha = .2, method = method, original_data = dat)
  coverage_df <- data.frame(
    lassoboot_ci,
    "truth" = dat$beta,
    "group" = j
  )

  coverages[[j]] <- coverage_df

}

coverages <- do.call(rbind, coverages)

save(coverages, file = glue("./rds/overall_v_feature_{dist}_{quantiles}_{method}_n{n}_p{p}.rds"))
