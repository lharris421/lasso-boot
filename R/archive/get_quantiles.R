#' Title
#'
#' @param quantile_functions holder
#' @param weights holder
#' @param len holder
#'
#' @return
#' @export
#'
#' @examples
get_quantiles <- function(quantile_functions, weights, len) {
  
  lengths <- round(len * weights)
  all_quantiles <- list()
  
  for (i in 1:length(weights)) {
    quant_func <- quantile_functions[[i]]
    all_quantiles[[i]] <- quant_func(1:lengths[i] / (lengths[i] + 1))
  }
  
  final_len <- sum(lengths)
  return(quantile(unlist(all_quantiles), 1:final_len / (final_len + 1)))
  
} 