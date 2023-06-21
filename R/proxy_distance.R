#' Title
#'
#' @param dist1 holder
#' @param dist2 holder
#'
#' @return
#' @export
#'
#' @examples
proxy_dist <- function(dist1, dist2) {
  mean(abs(sort(abs(dist1)) - sort(abs(dist2))))
}