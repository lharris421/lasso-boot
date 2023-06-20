#' Title
#'
#' @param dist1 
#' @param dist2 
#'
#' @return
#' @export
#'
#' @examples
rlt <- function(dist1, dist2) {
  move <- sort(abs(dist1)) - sort(abs(dist2))
  move_left <- abs(sum(move[move < 0]))
  move_right <- sum(move[move > 0])
  total_move <- sum(abs(move))
  return(list(move_left, move_right, total_move))
}
