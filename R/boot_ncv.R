#' Title
#'
#' @param X
#' @param y
#' @param cv_fit
#' @param lambda
#' @param sigma2
#' @param level
#' @param penalty
#' @param submethod
#' @param return_boot
#'
#' @return
#' @export
#'
#' @examples
boot_ncv <- function(X, y, cv_fit, lambda, sigma2,
                     enet_alpha = 1, gamma = switch(penalty, SCAD = 3.7, 3),
                     alpha = 0.05,
                     penalty = "lasso", submethod = "hybrid",
                     reselect_lambda = FALSE,
                     return_boot = FALSE, sigma2_reed = FALSE
                     ) {

  res <- boot_ncvreg(
    X, y, cv_fit, penalty = penalty, lambda = lambda, sigma2 = sigma2,
    alpha = enet_alpha, gamma = gamma,
    reselect_lambda = reselect_lambda, sigma2_reed = sigma2_reed
  )
  cis <- ci.boot_ncvreg(res, alpha = alpha, methods = submethod)

  cis <- cis %>%
    mutate(
      lambda = res$lambda,
      estimate = res$estimates,
      sigma = sqrt(res$sigma2)
    )

  if (return_boot) {
    list(cis = cis, boot = res)
  } else {
    cis
  }

}
