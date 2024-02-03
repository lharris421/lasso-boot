source("./scripts/setup/setup.R")

quantiles <- "disturbed"
method <- "quantile"
ns <- c(20, 30, 60)
p <- 30


plot_res <- list()
for (i in 1:length(ns)) {

  n <- ns[i]
  rt <- 2
  true_lambda <- (1 / n) * rt

  laplace_beta <- rlaplace(p, rate = rt)
  dat <- gen_data(n = n, p = p, beta = laplace_beta)
  dat$X <- ncvreg::std(dat$X)

  ### Selective Inference
  cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE)
  si_lam <- lam <- cv_res$lambda.min

  fit <- cv_res$glmnet.fit
  b <- coef(fit, s = lam, exact = TRUE)[-1]
  tryCatch({
    # sh <- estimateSigma(dat$X, dat$y)$sigmahat
    # print(sh)
    res <- fixedLassoInf(dat$X, dat$y, b, lam*length(dat$y), alpha = .2)
    B <- res$ci
    rownames(B) <- names(res$vars)
    colnames(B) <- c("lower", "upper")
    print(B)
    # B <- B[is.finite(B[,2]) & is.finite(B[,3]),-4]
    si_ci <- B %>%
      data.frame(method = "Selective Inference", variable = rownames(B)) %>%
      mutate(estimate = b[as.numeric(str_remove(rownames(B), "V"))])
  }, error = function(e) {
    si_ci <- NA
    print(e)
    }
  )

  ### HDI - Across a range of lambda values
  fit.lasso.allinfo <- boot.lasso.proj(dat$X, dat$y, return.bootdist = TRUE, B = nboot, boot.shortcut = TRUE)
  ci_hdi <- confint(fit.lasso.allinfo, level = 0.8)

  hdi_ci <- ci_hdi %>%
    data.frame(method = "BLP", variable = rownames(ci_hdi)) %>%
    mutate(estimate = fit.lasso.allinfo$bhat)
  hdi_lam <- fit.lasso.allinfo$lambda

  ### Lasso-boot
  lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = nboot, quantiles = quantiles)
  lassoboot_ci <- ci.boot.ncvreg(lassoboot, method = method, original_data = dat)
  lassoboot_lam <- lassoboot$lamdba

  plot_res[[i]] <- list(si_ci, hdi_ci, lassoboot_ci, si_lam, hdi_lam, lassoboot_lam, n, true_lambda)

}

save(plot_res, file = glue("./rds/method_comparison_laplace_{quantiles}_{method}.rds"))

