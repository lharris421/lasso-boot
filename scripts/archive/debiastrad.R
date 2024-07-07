source("./scripts/setup/setup.R")
library(tictoc)
library(mgcv)
library(patchwork)

## Function
# Find partials
calc_zj <- function(X, y, beta) {
  partial_residuals <-  y - (as.numeric(X %*% beta) - (X * matrix(beta, nrow=nrow(X), ncol=ncol(X), byrow=TRUE)))
  z <- (1/length(y))*colSums(X * partial_residuals)
  return(z)
}
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }
## Data arguments
corr <- "autoregressive"
rho <- 0.8
p <- 100
n <- p * 4
nboot <- 1000
simulations <- 100
alpha_orig <- .2
SNR <- 1

args_list <- list(method = "debiaslm",
                  data = "norm",
                  n = n,
                  snr = SNR,
                  correlation_structure = ifelse(rho > 0, corr, NA),
                  correlation = ifelse(!is.null(corr), rho * 100, NA),
                  nominal_coverage = alpha_orig * 100,
                  lambda = "cv",
                  p = p)

comb_res <- list()
coverages <- numeric(100)
averages <- numeric(100)
prop_lms <- numeric(100)
for (j in 1:100) {

  # dat <- hdrm::gen_data_abn()
  beta_vals <- rnorm(p)
  # beta_vals <- runif(p, -1, 1)
  dat <- gen_data_snr(n = n, p = p, p1 = p, beta = beta_vals, corr = corr, rho = rho, SNR = SNR)
  # beta <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
  # dat <- gen_data_snr(n = 100, p = 100, p1 = 100, beta = beta, rho = 0, SNR = 1)
  x_orig <- dat$X
  dat$X <- ncvreg::std(dat$X)
  rescale <- attr(dat$X, "scale")^(-1)

  dat$y <- dat$y - mean(dat$y)
  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso", alpha = 0.5)
  lambda <- cv_fit$lambda.min

  orig_betas <- coef(cv_fit$fit, lambda)[-1]
  orig_zjs <- calc_zj(dat$X, dat$y, orig_betas)

  # sigma2 <- cv_fit$cve[cv_fit$min] * (nrow(x_orig) / (nrow(x_orig) - sum(orig_betas!=0)))
  sigma2 <- cv_fit$cve[cv_fit$min]
  # sigma2 <- 1
  # sigma2 <- ((n - sum(lmin_beta[-1] != 0))^(-1)) * drop(crossprod(dat$y - (lmin_beta[1] + (dat$X %*% lmin_beta[-1]))))

  sd_norm_est <- numeric(ncol(dat$X))
  if (sum(orig_betas == 0) > 0) {
    sd_norm_est[orig_betas == 0] <- sqrt(sigma2 * (1 / diag(t(x_orig[,orig_betas == 0]) %*% x_orig[,orig_betas == 0])))
  }
  if (sum(orig_betas != 0) > 0) {
    sd_norm_est[orig_betas != 0] <- sqrt(sigma2 * diag(solve(t(x_orig[,orig_betas != 0]) %*% x_orig[,orig_betas != 0])))
  }

  zjs <- matrix(nrow = 1000, ncol = ncol(dat$X))
  betas <- matrix(nrow = 1000, ncol = ncol(dat$X))

  ## what about rel to average?
  prop_lm <- lambda / max(cv_fit$lambda)
  # alpha <- alpha_orig * (1 - mean(pmin(1, lambda / find_thresh(dat$X, dat$y)))) + ((alpha_orig / p) * mean(pmin(1, lambda / find_thresh(dat$X, dat$y))))

  print(prop_lm)

  prop_lms[j] <- prop_lm

  for (i in 1:nboot) {

    idx <- sample(1:length(dat$y), replace = TRUE)
    xnew <- dat$X[idx,]
    xnew <- ncvreg::std(xnew)
    rescaleNew <- attr(xnew, "scale")^(-1)
    ynew <- dat$y[idx]
    ynew <- ynew - mean(ynew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda <- lambda_max * prop_lm * 0.5
    lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
    if (lambda_min >= lambda_max | lambda >= lambda_max) { # Should review this
      lambda_max <- lambda + lambda / 100
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
    new_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)


    betanew <- coef(new_fit, lambda = lambda)[-1]
    betas[i,] <- betanew * (rescale * rescaleNew)

    cols <- which(betanew != 0)

    zj <- calc_zj(xnew, ynew, betanew)
    zjs[i,] <- zj




    # if (length(cols) > 1) {
    #   xsubset <- xnew[,cols]
    #   for (k in 1:length(cols)) {
    #     zjs[i,cols[k]] <- zjs[i,cols[k]] + sign(zjs[i,cols[k]])*lambda*summary(lm(xsubset[,k] ~ xsubset[,-k]))$r.squared
    #   }
    # }

    # if (length(cols) > 0) {
    #   zjs[i,cols] <- coef(lm(ynew ~ xnew[,cols]))[-1]
    # }

    # zjs[i,] <- zj * (rescale * rescaleNew)
    zjs[i,] <- zjs[i,] * rescaleNew


  }

  zjsnull <- matrix(nrow = 1000, ncol = ncol(dat$X))
  for (i in 1:nboot) {

    idx <- sample(1:length(dat$y), replace = FALSE)
    xnew <- dat$X
    xnew <- ncvreg::std(xnew)
    rescaleNew <- attr(xnew, "scale")^(-1)
    ynew <- dat$y[idx]
    ynew <- ynew - mean(ynew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda <- lambda_max * prop_lm
    lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
    if (lambda_min >= lambda_max | lambda >= lambda_max) { # Should review this
      lambda_max <- lambda + lambda / 100
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = 100))
    new_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    betanew <- coef(new_fit, lambda = lambda)[-1]

    zj <- calc_zj(xnew, ynew, betanew)
    zjsnull[i,] <- zj * rescaleNew


  }

  lowers <- numeric(ncol(zjs))
  uppers <- numeric(ncol(zjs))

  # norm_zjs <- abs(apply(zjs, 2, mean)) / sd_norm_est

  # prop_pen <- mean(pmin(1, lambda / find_thresh(dat$X, dat$y)))
  prop_pen <- mean(pmin(1, lambda / abs(orig_zjs)))

  pvals <- numeric(ncol(zjs))
  for (i in 1:ncol(zjs)) {
    pvals[i] <- mean(abs(zjsnull[,i]) > abs(orig_zjs[i]))
  }

  for (i in 1:ncol(zjs)) {

    # prop <- mean(betas[,i] == 0)
    # if (prop >= 0.8) {
    #
    #   curr_zjs <- zjs[betas[,i] == 0,i]
    #   lower_quantile <- (0.2 - (1 - (length(curr_zjs) / 1000))) / 2
    #   upper_quantile <- 1 - ((0.2 - (1 - (length(curr_zjs) / 1000))) / 2)
    #
    #   lower_quantile <- mean(zjs[,i] < quantile(curr_zjs, lower_quantile))
    #   uppper_quantile <- mean(zjs[,i] < quantile(curr_zjs, upper_quantile))
    #
    # } else if (prop < 0.8 & prop >= .2) {
    #
    #   lower_quantile <- 0.1
    #   upper_quantile <- 0.9
    #
    # } else if (prop < .2) {
    #
    #   adj <- (.2 - prop) * (0.5 - mean(zjs[betas[,i] != 0,i] < 0))
    #   lower_quantile <- 0.1 + adj
    #   upper_quantile <- 0.9 + adj
    #
    # }
    #
    # lower_quantile <- max(lower_quantile, (1/nrow(dat$X)))
    # upper_quantile <- min(upper_quantile, 1 - (1/nrow(dat$X)))
    # zjs_adj <- zjs[,i] * rescale[i]

    # lower_quantile <- (alpha / 2) - 0.1 * prop_lm[j]
    # upper_quantile <- (1 - (alpha / 2)) + 0.1 * prop_lm[j]

    # lower_quantile <- 0.1 - (0.1 * prop_lm)
    # upper_quantile <- 0.9 + (0.1 * prop_lm)


    # if (median(betas[,i]) == 0) {
    #   sample_points <- zjs[order(abs(zjs[,i])),i][1:800]
    #   left_points <- zjs[order(abs(zjs[,i])),i][801:1000]
    # } else if (median(betas[,i]) < 0) {
    #   sample_points <- zjs[order(zjs[,i]),i][1:800]
    #   left_points <- zjs[order(zjs[,i]),i][801:1000]
    # } else if (median(betas[,i]) > 0) {
    #   sample_points <- zjs[order(zjs[,i], decreasing = TRUE),i][1:800]
    #   left_points <- zjs[order(zjs[,i], decreasing = TRUE),i][801:1000]
    # } else {
    #   stop("stop")
    # }


    # sample_points <- sort(sample_points)
    #
    # ## Extend
    # if (min(sample_points) <= orig_betas[i] & max(sample_points) >= orig_betas[i]) {
    #   sample_points <- sample_points
    # } else if (orig_betas[i] < min(sample_points)) {
    #   # stop("stop")
    #   if (min(left_points) < orig_betas[i]) {
    #     new_lower <- max(left_points[left_points < orig_betas[i]])
    #   } else {
    #     new_lower <- min(left_points)
    #   }
    #
    #   to_add <- c(new_lower, left_points[left_points > new_lower])
    #
    #   sample_points <- c(to_add, sample_points[1:(800 - length(to_add))])
    #
    # } else if (orig_betas[i] > max(sample_points)) {
    #
    #   if (max(left_points) > orig_betas[i]) {
    #     new_upper <- min(left_points[left_points > orig_betas[i]])
    #   } else {
    #     new_upper <- max(left_points)
    #   }
    #
    #   to_add <- c(new_upper, left_points[left_points < new_upper])
    #
    #   sample_points <- c(to_add, sample_points[(length(to_add) + 1):800])
    #
    # }

    # sample_points <- updated_sample_points

    # if (length(sample_points) != 800) {stop("something wrong")}

    # if (sum(betas[,i] != 0) > 0) {
    #   shift <- 2 * lambda * mean(betas[,i] != 0) * (.5 - mean(betas[betas[,i] != 0,i] < 0))
    # } else {
    #   shift <- 0
    # }

    # curr_i <- which(order(orig_zjs, decreasing = TRUE)[1:100] == i)
    curr_i <- which(order(pvals)[1:100] == i)
    # prop_pen <- prop_lm
    alpha <- (alpha_orig * (1 - prop_pen)) + ((alpha_orig / (p - curr_i + 1)) * prop_pen)
    # alpha <- (alpha_orig * (1 - prop_pen)) + ((alpha_orig / p) * prop_pen)


    # nonzero_cert <- 1 - 2*pnorm(abs(mean(zjs[,i])) / sd_norm_est[i], lower.tail = FALSE)
    # nonzero_cert <- 1 - 2*pnorm(abs(mean(zjs[,i])) / sqrt(sigma2/n), lower.tail = FALSE)
    # nonzero_cert <- 1 - mean(abs(zjsnull[,i]) > abs(orig_zjs[i]))

    nonzero_cert <- 1 - pvals[i]
    # direction <- sign(mean(zjs[,i]))
    direction <- sign(orig_zjs[i])
    shift <- (alpha / 2) * nonzero_cert * direction

    lower_quantile <- (alpha / 2) + shift
    upper_quantile <- (1 - (alpha / 2))  + shift

    shift <- 0
    zjs_adj <- (zjs[,i] + shift) * rescale[i]
    # zjs_adj <- mean(zjs_adj) + ((zjs_adj - mean(zjs_adj)) * (sd_norm_est[i] / sd(zjs_adj)))

    lowers[i] <- quantile(zjs_adj, lower_quantile)
    uppers[i] <- quantile(zjs_adj, upper_quantile)

    # lowers[i] <- min(sample_points) * rescale[i]
    # uppers[i] <- max(sample_points) * rescale[i]


  }

  res <- data.frame(truth = dat$beta, lower = lowers, upper = uppers, group = j)

  comb_res[[j]] <- res

  coverages[j] <- res %>%
    mutate(covered = lower <= truth & upper >= truth) %>%
    pull(covered) %>%
    mean()

  print(glue::glue("Actual coverage: {coverages[j]}"))
  print(glue::glue("Cumulative coverage: {mean(coverages[1:j])}"))


}

per_var_data <- do.call(rbind, comb_res)

per_var_data %>%
  mutate(absb = abs(truth), covered = lower <= truth & upper >= truth) %>%
  group_by(absb) %>%
  dplyr::summarise(coverage = mean(covered))

res_list <- list("per_var_data" = per_var_data)
## Need to update to take grid or way to subset list easily
save_objects(folder = rds_path, list("per_var_data" = per_var_data), args_list = args_list, overwrite = TRUE, save_method = "rds")

plots <- list()

cutoff <- 0.275
ns <- unique(per_var_data$n)

# Function to calculate model results
calculate_model_results <- function(data) {
  data %>%
    mutate(
      covered = lower <= truth & upper >= truth,
      mag_truth = abs(truth),
      covered = as.numeric(covered)
    )
}

# Function to perform fitting and prediction
predict_covered <- function(data, x_values, method) {
  fit <- gam(covered ~ s(mag_truth) + s(group, bs = "re"), data = data, family = binomial)
  y_values <- predict(fit, data.frame(mag_truth = x_values, group = 101), type = "response")
  data.frame(x = x_values, y = y_values, method = method)
}

model_res <- calculate_model_results(per_var_data)

xvals <- seq(from = 0, to = cutoff, length.out = cutoff * 100 + 1)
density_data <- data.frame(x = xvals, density = 2 * dlaplace(xvals, rate = 14.14))

tmp <- model_res

cat("Average coverage: ", mean(tmp$covered), "\n")

xs <- seq(0, cutoff, by = 0.01)
line_data <- predict_covered(tmp, xs, "debiased")
line_data_avg <- data.frame(avg = mean(tmp$covered), method = "debiased")

ggplot() +
  geom_line(data = line_data %>% mutate(method = "debiased"), aes(x = x, y = y, color = method)) +
  geom_hline(data = line_data_avg, aes(yintercept = avg, color = method), linetype = 2) +
  geom_hline(aes(yintercept = 0.8), linetype = 1, alpha = .5) +
  theme_bw() +
  xlab(expression(abs(beta))) +
  ylab("Coverage") +
  coord_cartesian(ylim = c(0, 1))

## lhb
miss_low <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    sign(truth) == 0 ~ 0,
    sign(truth) == 1 & sign(upper) == 1 & upper < truth  ~ 1,
    sign(truth) == -1 & sign(lower) == -1 & lower > truth ~ 1,
    sign(truth) != 0 & sign(lower) == 0 & sign(upper) == 0 ~ 1,
    TRUE ~ 0
  )
}
miss_high <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    sign(truth) == 1 & lower > truth ~ 1,
    sign(truth) == -1 & upper < truth ~ 1,
    TRUE ~ 0
  )
}
sign_inversion <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    (sign(truth) == 1 & sign(upper) == -1) | (sign(truth) == 1 & sign(upper) == 0 & sign(lower) == -1) ~ 1,
    (sign(truth) == -1 & sign(lower) == 1) | (sign(truth) == -1 & sign(lower) == 0 & sign(upper) == 1) ~ 1,
    TRUE ~ 0
  )
}

xs <- seq(0, cutoff, by = .01)
plot_res <- list()

plot_data <- per_var_data %>%
  rowwise() %>%
  mutate(bias_high = miss_high(lower, upper, truth)) %>%
  mutate(bias_low = miss_low(lower, upper, truth)) %>%
  mutate(bias_sign = sign_inversion(lower, upper, truth)) %>%
  mutate(mag_truth = abs(truth))

fit <- gam(bias_high ~ s(mag_truth) + s(group, bs = "re"), data = plot_data, family = binomial)
ys_high <- predict(fit, data.frame(mag_truth = xs, group = 101), type ="response")
fit <- gam(bias_low ~ s(mag_truth) + s(group, bs = "re"), data = plot_data, family = binomial)
ys_low <- predict(fit, data.frame(mag_truth = xs, group = 101), type ="response")
fit <- gam(bias_sign ~ s(mag_truth) + s(group, bs = "re"), data = plot_data, family = binomial)
ys_sign <- predict(fit, data.frame(mag_truth = xs, group = 101), type ="response")

plot_res <- data.frame(xs = xs, bias = ys_low - ys_high, bias_low = ys_low, bias_high = ys_high, bias_sign = ys_sign)

plot_bias_h <- plot_res %>%
  ggplot(aes(x = xs, y = bias_high)) +
  geom_line() +
  theme_bw() +
  xlab(NULL) +
  ylab("P(Miss Away)") +
  theme(legend.position = "none")

plot_bias_l <- plot_res %>%
  ggplot(aes(x = xs, y = bias_low + bias_sign)) +
  geom_line() +
  theme_bw() +
  xlab(NULL) +
  ylab("P(Miss Towards + T3)") +
  theme(legend.position = "none")


# Create a layout matrix
(plot_bias_l + plot_bias_h) +
  plot_layout(guides = "collect", axis_titles = "collect")

