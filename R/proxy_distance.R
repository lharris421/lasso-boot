#' Title
#'
#' @param dist1 
#' @param dist2 
#'
#' @return
#' @export
#'
#' @examples
proxy_dist <- function(dist1, dist2) {
  mean(abs(sort(abs(dist1)) - sort(abs(dist2))))
}