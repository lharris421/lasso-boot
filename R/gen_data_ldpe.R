#' Title
#'
#' @param n
#' @param p
#' @param a
#' @param rho
#' @param corr
#' @param which_nonzero
#'
#' @return
#' @export
#'
#' @examples
gen_data_ldpe <- function(n, p, a, rho, corr = "autoregressive",
                          which_nonzero = c(1500, 1800, 2100, 2400, 2700, 3000)) {

  true_lambda <- sqrt((2/n)*log(p))
  betas <- 3*true_lambda / ((1:p)^a)
  betas[which_nonzero] <- 3*true_lambda

  gen_data(n, p, beta = betas, corr = corr, rho = rho)

}
