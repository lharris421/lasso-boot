#' Title
#'
#' @param X
#' @param y
#' @param return.bootdist
#' @param B
#' @param boot.shortcut
#' @param alpha
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
blp <- function(X, y, return.bootdist = TRUE, B = 1000, boot.shortcut = FALSE, alpha = 0.2, lambda = NULL) {

  tryCatch({
    fit.lasso.allinfo <- boot.lasso.proj(
      X, y,
      return.bootdist = return.bootdist,
      lambda = lambda,
      B = B,
      boot.shortcut = boot.shortcut,
      verbose = FALSE
    )
    ci_hdi <- confint(fit.lasso.allinfo, level = 1 - alpha)

    ci <- ci_hdi %>%
      data.frame(variable = rownames(ci_hdi)) %>%
      mutate(estimate = fit.lasso.allinfo$bhat) %>%
      dplyr::select(variable, estimate, lower, upper) %>%
      mutate(
        lambda = fit.lasso.allinfo$lambda,
        sigma = fit.lasso.allinfo$sigmahat
      )

    res <- ci
    return(res)
  }, error = function(e) {
    print(e);
    res <- data.frame(variable = colnames(X), estimate = NA, lower = NA, upper = NA, lambda = NA, sigma = NA, error = as.character(e))
    return(res)
  })
}
