blp <- function(dat, return.bootdist = TRUE, B = 1000, boot.shortcut = FALSE, alpha = 0.2, lambda = NULL) {

  tryCatch({
    fit.lasso.allinfo <- boot.lasso.proj(
      dat$X, dat$y,
      return.bootdist = return.bootdist,
      lambda = lambda,
      B = B,
      boot.shortcut = boot.shortcut,
      verbose = TRUE
    )
    ci_hdi <- confint(fit.lasso.allinfo, level = 1 - alpha)

    ci <- ci_hdi %>%
      data.frame(method = "blp", variable = rownames(ci_hdi)) %>%
      mutate(estimate = fit.lasso.allinfo$bhat)
    lam <- fit.lasso.allinfo$lambda
    sig <- fit.lasso.allinfo$sigma
    res <- list("confidence_interval" = ci, "lambda" = lam, "sigma" = sig)
    return(res)
  }, error = function(e) {
    print(e);
    res <- NULL
    return(res)
  })
}
