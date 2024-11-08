#' Title
#'
#' @param method_info
#' @param data
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
run_data_analysis <- function(method_info, data, overwrite = TRUE) {

  # Combine the method information and simulation information
  analysis_arguments <- c(
    method_info, data
  )

  # Run the simulation
  analysis_results <- do.call(data_analysis, analysis_arguments)
  analysis_results$call[1] <- "data_analysis" ## Need more elegant solution

  # Save the results
  indexr::save_objects("./rds", analysis_results, overwrite = overwrite)

}
data_analysis <- function(method, method_arguments = list(), alpha = 0.05, data) {

  original_call <- match.call()

  results <- do.call(method, c(list(X = data$X, y = data$y, alpha = alpha), method_arguments))

  return(list("call" = original_call, "results" = results))

}
