true_rate <- sqrt(2*p)
true_lambda <- true_rate / n

lambda_max <- max(ncvreg(dat$X, dat$y, penalty = method)$lambda)
lambda_min <- (true_lambda / lambda_max) * .999
if (modifier %in% c("tl", "tls")) {
  lambda <- true_lambda
  if (modifier == "tl") {
    ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
  } else if (modifier == "tls") {
    sigma2 <- 1
  }
}

cv_fit <- cv.ncvreg(
  dat$X, dat$y, penalty = penalty,
  lambda.min = ifelse(!(modifier %in% c("tl", "tls")), ifelse(n > p, 0.001, 0.05), lambda_min),
  max.iter = 1e8, alpha = enet_alpha, gamma = gamma
)


lassoboot <- boot_ncvreg(
  dat$X, dat$y, penalty = method, verbose = FALSE, nboot = nboot,
  max.iter = 1e8, lambda = lambda, sigma2 = sigma2, lambda.min = lambda_min,
  alpha = enet_alpha, gamma = gamma)
ci <- ci.boot_ncvreg(lassoboot, alpha = alpha)
