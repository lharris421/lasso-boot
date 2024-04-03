plot_ridge <- function(ridge_ci, n = 20, quiet = TRUE) {

  plot_res <- ridge_ci[-1,] %>%
    data.frame()

  plot_res$variable <- rownames(plot_res)
  plot_res <- plot_res %>%
    dplyr::arrange(desc(abs(estimate))) %>%
    head(n)

  plot_res$variable <- factor(plot_res$variable, levels = rev(plot_res$variable))

  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable)) +
    geom_point(aes(x = estimate, y = variable)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")

}
