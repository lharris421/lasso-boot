library(glmnet)
library(tictoc)
library(selectiveInference)
library(hdi)
.libPaths("./local")
library(ncvreg)
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

res_coverage <- res_time <- res_lambda <- matrix(nrow = 3, ncol = 3)
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
    tryCatch({
      si_coverage <- logical(p)
      start <- Sys.time()
      cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE)
      lam <- cv_res$lambda.min

      fit <- cv_res$glmnet.fit
      b <- coef(fit, s = lam, exact = TRUE)[-1]
      sh <- estimateSigma(dat$X, dat$y)$sigmahat
      res <- fixedLassoInf(dat$X, dat$y, b, lam*length(dat$y), sigma=sh, alpha = .2)
      times[j,1] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      bb <- res$vmat %*% dat$y
      B <- cbind(bb, res$ci, res$pv)
      rownames(B) <- names(res$vars)
      # B <- B[is.finite(B[,2]) & is.finite(B[,3]),-4]
      si_ci <- B %>%
        data.frame(method = "Selective Inference", variable = rownames(B)) %>%
        rename(estimate = X1, lower = X2, upper = X3) %>%
        left_join(truth_df, by = "variable")

      si_coverage[as.numeric(str_remove(si_ci$variable, "V"))] <- si_ci$lower <= si_ci$truth & si_ci$truth <= si_ci$upper
      lambdas[j, 1] <- lam
    }, error = function(e) {print(e); print("SI failed"); si_ci <- NULL; lambdas[j, 1] <- NA})

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

    ### Lasso-boot
    start <- Sys.time()
    lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE)
    times[j,3] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    lassoboot_ci <- ci.boot.ncvreg(lassoboot)
    lassoboot_coverage <- lassoboot_ci$lower <= laplace_beta & laplace_beta <= lassoboot_ci$upper
    lambdas[j, 3] <- lassoboot$lamdba

    ## Coverage calculation
    coverage_df <- cbind(si_coverage, hdi_coverage, lassoboot_coverage)
    coverages[[j]] <- coverage_df
    overall_coverage[j,] <- apply(coverage_df, 2, mean)

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

  all_coverages[[i]] <- coverages

}

