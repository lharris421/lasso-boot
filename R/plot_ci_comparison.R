plot_ci_comparison <- function(cis, nvars = 20) {

  plot_vars <- list()
  for (current_method in unique(cis$method)) {
    plot_vars[[current_method]] <- cis %>%
      filter(method == current_method) %>%
      dplyr::arrange(desc(abs(estimate))) %>%
      slice_head(n = nvars) %>%
      dplyr::arrange(desc(estimate)) %>%
      pull(variable)
  }
  plot_vars <- unique(unlist(plot_vars))

  plot_res <- cis %>%
    filter(variable %in% plot_vars) %>%
    dplyr::arrange(desc(estimate))

  plot_res$variable <- factor(plot_res$variable, levels = rev(plot_vars))

  plot_res %>%
    mutate(method = factor(methods_pretty[method], levels = methods_pretty[unique(cis$method)])) %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable)) +
    geom_point(aes(x = estimate, y = variable)) +
    theme_bw() +
    ylab(NULL) + xlab(NULL) +
    scale_color_manual(name = "Method", values = colors) +
    theme(legend.position = "none",
          # legend.position = c(.9, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.background = element_rect(fill = "transparent")) +
    facet_wrap(~method, scales = "free_x", nrow = 1)

}
