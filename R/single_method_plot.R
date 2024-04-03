single_method_plot <- function(per_var_data, ns, alpha) {
  cutoff <- round(quantile(abs(per_var_data$truth), .98), 1)
  line_data <- list()
  line_data_avg <- list()
  for (j in 1:length(ns)) {

    tmp <- per_var_data %>%
      mutate(
        covered = lower <= truth & upper >= truth,
        mag_truth = abs(truth), covered = as.numeric(covered)
      ) %>%
      filter(n == ns[j])

    fit <- gam(covered ~ s(mag_truth) + s(group, bs = "re"), data = tmp, family = binomial)
    xs <- seq(0, cutoff, by = .01)
    ys <- predict(fit, data.frame(mag_truth = xs, group = 101), type ="response")

    line_data[[j]] <- data.frame(x = xs, y = ys, n = factor(ns[j]))
    line_data_avg[[j]] <- data.frame(avg = mean(tmp$covered), n = factor(ns[j]))
  }

  line_data_avg <- do.call(rbind, line_data_avg)
  line_data <- do.call(rbind, line_data)
  nom_data <- data.frame(alpha = alpha)

  plt <- ggplot() +
    geom_line(data = line_data, aes(x = x, y = y, color = n)) +
    geom_hline(data = line_data_avg, aes(yintercept = avg, color = n), linetype = 2) +
    geom_hline(data = nom_data, aes(yintercept = (1 - alpha)), linetype = 1, alpha = .5) +
    geom_text() +
    theme_bw() +
    xlab(expression(abs(beta))) +
    ylab(NULL) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, round(quantile(abs(per_var_data$truth), .98), 1))) +
    scale_color_manual(name = "N", values = colors)

  return(plt)
}
