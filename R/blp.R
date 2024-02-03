blp <- function(dat, alpha = 0.2) {

  tryCatch({
    fit.lasso.allinfo <- boot.lasso.proj(dat$X, dat$y, return.bootdist = TRUE, B = nboot, boot.shortcut = TRUE)
    ci_hdi <- confint(fit.lasso.allinfo, level = 1 - alpha)

    ci <- ci_hdi %>%
      data.frame(method = "blp", variable = rownames(ci_hdi)) %>%
      mutate(estimate = fit.lasso.allinfo$bhat)
    lam <- fit.lasso.allinfo$lambda
    sig <- fit.lasso.allinfo$sigma
    res <- list("confidence_interval" = ci, "lambda" = lam, "sigma" = sigma)
    return(res)
  }, error = function(e) {
    print(e);
    res <- NULL
    return(res)
  })
}
