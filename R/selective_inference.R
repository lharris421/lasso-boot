selective_inference <- function(dat, alpha = .2, estimate_sigma = FALSE, sigma = NULL) {
  res <- list()
  # cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE, thresh = 1e-10)
  cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE)
  lam <- cv_res$lambda.min
  fit <- cv_res$glmnet.fit
  # b <- coef(fit, x=dat$X, y=dat$y, s = lam, exact = TRUE)[-1]
  b <- coef(fit, s = lam)[-1]
  tryCatch({
    if (estimate_sigma & is.null(sigma)) {
      sigma <- estimateSigma(dat$X, dat$y)$sigmahat
    }
    res <- fixedLassoInf(dat$X, dat$y, b, lam*length(dat$y), alpha = alpha, sigma = sigma)
    B <- res$ci
    rownames(B) <- names(res$vars)
    colnames(B) <- c("lower", "upper")
    # B <- B[is.finite(B[,2]) & is.finite(B[,3]),-4]
    estimates <- data.frame(variable = colnames(dat$X), estimate = b)
    scale_df <- data.frame(variable = colnames(dat$X), scale = attr(dat$X, "scale"))
    si_ci <- B %>%
      data.frame(method = "selective_inference", variable = rownames(B)) %>%
      left_join(estimates) %>%
      left_join(scale_df) %>%
      mutate(estimate = estimate / scale, lower = lower / scale, upper = upper / scale)
    res <- list("confidence_interval" = si_ci, "lambda" = lam, "sigma" = sigma)
    return(res)
  }, error = function(e) {
    print(e);
    res <- NULL
    return(res)
  }
  )
}
