#' Title
#'
#' @param lam 
#' @param n 
#' @param sig2 
#'
#' @return
#' @export
#'
#' @examples
rt <- function(lam, n, sig2 = 1) {
  (n*lam) / sig2
}
