#' Title
#'
#' @param pzero 
#' @param rate 
#' @param n 
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