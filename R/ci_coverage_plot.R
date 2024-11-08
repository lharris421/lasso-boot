#' Title
#'
#' @param plot_list
#' @param vars
#'
#' @return
#' @export
#'
#' @examples
ci_coverage_plot <- function(plot_list, vars) {

  results <- do.call(bind_rows, plot_list) %>%
    filter(variable %in% vars) %>%
    mutate(covered = lower <= truth & truth <= upper)

  coverage_labels <- results %>%
    group_by(method, variable) %>%
    summarise(coverage = mean(covered)) %>%
    mutate(label = round(coverage, 3))

  results %>%
    ggplot() +
    geom_errorbar(aes(ymin = lower, ymax = upper, x = iteration, color = covered), alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) +
    geom_hline(aes(yintercept = truth), color = "gold", linetype = "dashed") +
    theme_minimal() +
    xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_text(
      data = coverage_labels,
      aes(x = -Inf, y = Inf, label = label),
      hjust = 0, vjust = 1,
      inherit.aes = FALSE,
      size = 2
    ) +
    facet_grid(method ~ variable)

}
