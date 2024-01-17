rm(list=ls())
unloadNamespace("ncvreg")
.libPaths("./local")
library(ncvreg)
library(glmnet)
library(tictoc)
library(selectiveInference)
library(hdi)
library(ggplot2)
library(hdrm)
library(dplyr)
library(stringr)

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}

my_seed <- 189807771
set.seed(my_seed)

ns <- c(20, 30, 60)
p <- 30
rt <- 2

res_width <- res_coverage <- res_time <- res_lambda <- matrix(nrow = 3, ncol = 3)
all_coverages <- list()

for (i in 1:length(ns)) {

  print(i)
  lambdas <- matrix(nrow = 100, ncol = 3)
  times <- matrix(nrow = 100, ncol = 3)
  overall_coverage <- matrix(nrow = 100, ncol = 3) ## make back to 3
  widths <- matrix(nrow = 100, ncol = 3)
  coverages <- list()

  n <- ns[i]
  true_lambda <- (1 / n) * rt

  for (j in 1:100) {

    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data(n = n, p = p, beta = laplace_beta)
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    ### Lasso-boot - original
    start <- Sys.time()
    lassoboot_o <- boot.ncvreg(dat$X, dat$y, verbose = FALSE)
    times[j,1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    lassoboot_ci_o <- ci.boot.ncvreg(lassoboot_o)
    lassoboot_coverage_o <- lassoboot_ci_o$lower <= laplace_beta & laplace_beta <= lassoboot_ci_o$upper
    lambdas[j,1] <- lassoboot_o$lamdba
    widths[j,1] <- mean(lassoboot_ci_o$upper - lassoboot_ci_o$lower)

    ### Lasso-boot - combine
    start <- Sys.time()
    lassoboot_c <- boot.ncvreg.r(dat$X, dat$y, verbose = FALSE, quantiles = c(0.1, 0.9))
    times[j,2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    lassoboot_ci_c <- ci.boot.ncvreg.r(lassoboot_c)
    lassoboot_coverage_c <- lassoboot_ci_c$lower <= laplace_beta & laplace_beta <= lassoboot_ci_c$upper
    lambdas[j, 2] <- lassoboot_c$lamdba
    widths[j,2] <- mean(lassoboot_ci_c$upper - lassoboot_ci_c$lower)

    ### Lasso-boot - sample
    start <- Sys.time()
    lassoboot_s <- boot.ncvreg(dat$X, dat$y, verbose = FALSE)
    times[j,3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    lassoboot_ci_s <- ci.boot.ncvreg(lassoboot_s)
    lassoboot_coverage_s <- lassoboot_ci_s$lower <= laplace_beta & laplace_beta <= lassoboot_ci_s$upper
    lambdas[j, 3] <- lassoboot_s$lamdba
    widths[j,3] <- mean(lassoboot_ci_s$upper - lassoboot_ci_s$lower)

    ## Coverage calculation
    # coverage_df <- cbind(si_coverage, hdi_coverage, lassoboot_coverage)
    coverage_df <- data.frame(
      "lasso_original" = lassoboot_coverage_o,
      "lasso_combined" = lassoboot_coverage_c,
      "lasso_sample" = lassoboot_coverage_s,
      "truth" = laplace_beta, "group" = j
    )

    coverages[[j]] <- coverage_df
    overall_coverage[j,] <- apply(coverage_df[,-c((ncol(coverage_df)-1) : ncol(coverage_df)),drop=FALSE], 2, mean)

  }

  mean_times <- apply(times, 2, mean, na.rm = TRUE)
  sd_times <- apply(times, 2, sd, na.rm = TRUE)
  res_time[i,] <- paste0(round(mean_times, 3), " (", round(sd_times, 3), ")")

  mean_coverages <- apply(overall_coverage, 2, mean)
  sd_coverages <- apply(overall_coverage, 2, sd)
  res_coverage[i,] <- paste0(round(mean_coverages, 3), " (", round(sd_coverages, 3), ")")

  mean_lambdas <- apply(lambdas, 2, mean, na.rm = TRUE)
  sd_lambdas <- apply(lambdas, 2, sd, na.rm = TRUE)
  res_lambda[i,] <- paste0(round(mean_lambdas, 3), " (", round(sd_lambdas, 3), ")")

  mean_widths <- apply(widths, 2, mean, na.rm = TRUE)
  sd_widths <- apply(widths, 2, sd, na.rm = TRUE)
  res_width[i,] <- paste0(round(mean_widths, 3), " (", round(sd_widths, 3), ")")

  all_coverages[[i]] <- coverages

}

save(res_time, res_coverage, res_lambda, res_width, all_coverages, file = "./rds/method_comparison_laplace_100.rds")
