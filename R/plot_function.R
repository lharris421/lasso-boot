plot_function <- function(plot_list) {

  lambda_max <- 1
  lambda_min <- 0.001
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 10))

  cv_coverage <- plot_list %>%
    mutate(covered = truth >= lower & truth <= upper)  %>%
    filter(lambda == 11) %>%
    pull(covered) %>%
    mean()

  plot_list <- plot_list %>%
    mutate(covered = truth >= lower & truth <= upper)  %>%
    filter(lambda <= 10)

  overall_cov <- plot_list %>%
    dplyr::mutate(lambda = lambda_seq[lambda]) %>%
    dplyr::group_by(lambda) %>%
    dplyr::summarise(coverage = mean(covered)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = "Overall")
  size_cov <- plot_list %>%
    dplyr::mutate(lambda = lambda_seq[lambda]) %>%
    dplyr::mutate(group = as.character(group)) %>%
    dplyr::group_by(lambda, group) %>%
    dplyr::summarise(coverage = mean(covered)) %>%
    dplyr::ungroup()
  coverage_data <- dplyr::bind_rows(overall_cov, size_cov) %>%
    dplyr::mutate(group = factor(group, levels = c("Small", "Moderate", "Large", "Overall")))

  gg <- ggplot(data = coverage_data, aes(x = lambda, y = coverage, group = group, color = group)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linewidth = .5, linetype = 2) +
    theme_bw() +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks('log10', function(x) 10^x),
                       labels = trans_format('log10', math_format(10^.x))) +
    coord_cartesian(xlim = c(1, .001), ylim = c(0, 1.0)) +
    scale_color_manual(name = expression(abs(beta)), values = colors) +
    ggtitle(glue("{plot_list$dist_type}"), subtitle = glue("CV Coverage: {cv_coverage}"))

  return(gg)

}
