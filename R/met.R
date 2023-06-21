#' Title
#'
#' @param lam holder
#' @param sig2 holder
#' @param n holder
#' @param truth holder
#'
#' @return
#' @export
#'
#' @examples
met <- function(lam, sig2 = 1, n, truth) {
  rate <- (n*lam) / sig2
  prior <- qexp((1:length(truth)) / (length(truth) + 1), rate = rate)
  return(proxy_dist(prior, truth))
}