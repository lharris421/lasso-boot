rm(list = ls())
source("./scripts/setup.R")
library(ncvreg)

data <- hdrm::gen_data(family = "poisson", n = 200, p = 50, p1 = 10)
data$X <- ncvreg::std(data$X)
pipe(data$X, data$y, family = "poisson", penalty = "lasso") %>%
  mutate(
    truth = data$beta,
    covered = lower <= truth & truth <= upper,
    abs_diff = abs(coef - estimate),
    lambda_realized = lambda / weights
  )

# beta_seq <- seq(-0.5, 0.5, by = .001)
# lambda <- 0.08915687
# se <- 0.1101055
# zjs <- c(0.26250709)
#
# post_dens <- function(x, lambda, zj, se) {
#   if (x < 0) {
#     exp((zj*lambda)/se)*exp(-(1 / (2*se))*(x - (zj + lambda))^2)
#   } else {
#     exp(-(zj*lambda)/se)*exp(-(1 / (2*se))*(x - (zj - lambda))^2)
#   }
# }
#
# res <- list()
# for (i in 1:length(zjs)) {
#   dens <- sapply(beta_seq, post_dens, lambda, zjs[i], se)
#   # dens <- dens / max(dens)
#   res[[i]] <- data.frame(x = beta_seq, y = dens, zj = zjs[i])
# }
# dat <- do.call(dplyr::bind_rows, res)
#
#
# ggplot2::ggplot() +
#   ggplot2::geom_area(data = dat, aes(x = x, y = y), alpha = 0.5) +
#   ggplot2::facet_wrap(~zj, labeller = label_bquote(z == .(zj))) +
#   theme_bw() +
#   theme(panel.spacing = unit(1, "lines")) +
#   xlab(expression(beta)) +
#   ylab("Density") +
#   geom_vline(xintercept = 0.1733502165) +
#   geom_vline(xintercept = beta_seq[which.max(dens)], color = "red")
