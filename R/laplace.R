#' Title
#'
#' @param n
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}
#' Title
#'
#' @param x
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
dlaplace <- function(x, rate = 1) {
  dexp(abs(x), rate) / 2
}
