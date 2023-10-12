#' Title
#'
#' @param lam holder
#' @param n holder
#' @param sig2 holder
#'
#' @return
#' @export
#'
#' @examples
rt <- function(lam, n, sig2 = 1) {
  (n*lam) / sig2
}
