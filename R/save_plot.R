#' Title
#'
#' @param plot
#' @param file_path
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
save_plot <- function(plot, file_path, width = 6, height = 4) {
  # Open a PDF device to save the plot
  pdf(file_path, width = width, height = height)

  # Print the plot to the device
  print(plot)

  # Close the device
  dev.off()

  # Return the path to the saved file
  return(file_path)
}
