#' Title
#'
#' @param res_list
#' @param truth
#' @param vars
#'
#' @return
#' @export
#'
#' @examples
compare_methods <- function(res_list, truth = NULL, vars = NULL) {

  # Combine all results into one data frame, using the list names as method names
  combined_res <- purrr::imap_dfr(res_list, ~ .x %>%
                                    dplyr::select(variable, estimate, lower, upper) %>%
                                    dplyr::mutate(method = .y))

  if (!is.null(vars)) {
    combined_res <- combined_res %>%
      filter(variable %in% vars)
  }

  # Create the plot
  plot <- ggplot(combined_res) +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable, color = method), alpha = 0.8) +
    geom_point(aes(x = estimate, y = variable, color = method), alpha = 0.8) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")

  if (!is.null(truth)) {

    if (!is.null(vars)) {
      truth <- truth %>%
        filter(var %in% vars)
    }

    plot <- plot +
      geom_point(data = truth, aes(y = var, x = truth), color = "red", shape = 1)
  }

  return(plot)

}

