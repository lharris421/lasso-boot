#' Title
#'
#' @param pzero holder
#' @param rate holder
#' @param n holder
#'
#' @return
#' @export
#'
#' @examples
discrete_mixture <- function(pzero = .8, rate = 1, n = 1e5) {
  
  zeros <- rep(0, round(n*pzero))
  nonzero <- n - length(zeros)
  exp_quantiles <- qexp(1:nonzero / (nonzero + 1), rate = rate)
  return(c(zeros, exp_quantiles))
  
}