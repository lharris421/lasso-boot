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
  overall_coverage <- matrix(nrow = 100, ncol = 3)
  widths <- matrix(nrow = 100, ncol = 3) ## Need to ADD
  coverages <- list()

  n <- ns[i]
  true_lambda <- (1 / n) * rt

  for (j in 1:100) {
    laplace_beta <- rlaplace(p, rate = rt)
    dat <- gen_data(n = n, p = p, beta = laplace_beta)
    dat$X <- ncvreg::std(dat$X)
    truth_df <- data.frame(variable = names(dat$beta), truth = dat$beta)

    ### Selective Inference
    si_coverage <- logical(p)
    start <- Sys.time()
    cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE)
    lam <- cv_res$lambda.min
    fit <- cv_res$glmnet.fit
    b <- coef(fit, x=dat$X, y=dat$y, s = lam, exact = TRUE)[-1]
    tryCatch({
      # sh <- estimateSigma(dat$X, dat$y)$sigmahat
      res <- fixedLassoInf(dat$X, dat$y, b, lam*length(dat$y), alpha = .2)
      times[j,1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      B <- res$ci
      rownames(B) <- names(res$vars)
      colnames(B) <- c("lower", "upper")
      # B <- B[is.finite(B[,2]) & is.finite(B[,3]),-4]
      si_ci <- B %>%
        data.frame(method = "Selective Inference", variable = rownames(B)) %>%
        mutate(estimate = b[as.numeric(str_remove(rownames(B), "V"))]) %>%
        left_join(truth_df, by = "variable")

      si_coverage[as.numeric(str_remove(si_ci$variable, "V"))] <- si_ci$lower <= si_ci$truth & si_ci$truth <= si_ci$upper
      lambdas[j, 1] <- lam
      widths[j,1] <- mean(si_ci$upper - si_ci$lower) ## Only for variables included
    }, error = function(e) {
      print(e); print("SI failed")
      si_ci <- NULL
      }
    )

    ### HDI - Across a range of lambda values
    start <- Sys.time()
    fit.lasso.allinfo <- boot.lasso.proj(dat$X, dat$y, return.bootdist = TRUE, B = 100, boot.shortcut = TRUE)
    times[j,2] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    ci_hdi <- confint(fit.lasso.allinfo, level = 0.8)

    hdi_ci <- ci_hdi %>%
      data.frame(method = "BLP", variable = rownames(ci_hdi)) %>%
      mutate(estimate = fit.lasso.allinfo$bhat)
    hdi_coverage <- hdi_ci$lower <= laplace_beta & laplace_beta <= hdi_ci$upper
    lambdas[j, 2] <- fit.lasso.allinfo$lambda
    widths[j,2] <- mean(hdi_ci$upper - hdi_ci$lower)

    ### Lasso-boot
    start <- Sys.time()
    lassoboot <- boot.ncvreg.r(dat$X, dat$y, verbose = FALSE)
    times[j,3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    lassoboot_ci <- ci.boot.ncvreg.r(lassoboot)
    lassoboot_coverage <- lassoboot_ci$lower <= laplace_beta & laplace_beta <= lassoboot_ci$upper
    lambdas[j, 3] <- lassoboot$lamdba
    widths[j,3] <- mean(lassoboot_ci$upper - lassoboot_ci$lower)

    coverage_df <- data.frame(
      "si" = si_coverage,
      "hdi" = hdi_coverage,
      "lasso_sample" = lassoboot_coverage,
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
