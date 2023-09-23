## Laplace functions
dlaplace <- function(x, rate = 1) {
  dexp(abs(x), rate) / 2
}

qlaplace <- function(p, rate = 1) {
  if (p <= .5) {
    p <- (.5 - p)*2
    qexp(p, rate)*-1
  } else {
    p <- (p - .5)*2
    qexp(p, rate)
  }
}

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}

find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }

## Log-liklihood using partial residuals
ll <- function(beta, partial_residuals, sigma, xvar, type = "original", ynew) {

  # Compute the log-likelihood(s)
  if (type == "univariate") {
    ## z <- (1/length(partial_residuals))*as.numeric(t(xvar) %*% partial_residuals)
    z <- (1/length(partial_residuals))*as.numeric(t(xvar) %*% ynew)
    tmp <- sapply(beta, function(x) dnorm(z - x, mean = 0, sd = sigma, log = TRUE))
  } else {
    tmp <- sapply(beta, function(x) sum(dnorm(partial_residuals - xvar*x, mean = 0, sd = sigma, log = TRUE)))
  }

  return(tmp)

}

## Computing density for a given beta for the posterior
density_function <- function(x, rate, partial_residuals, sigma, xvar, normalizer = 1, multiplier = 0, type = "original", log = FALSE, ynew) {

  prior <- log(dlaplace(x, rate = rate))
  llik <- ll(beta = x, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar, type = type, ynew = ynew)

  if (!log) {
    ret <- exp(prior + llik - log(normalizer) + multiplier)
  } else {
    ret <- prior + llik - log(normalizer) + multiplier
  }
  return(ret)

}

## Objective
obj <- function(beta, p, sigma, rate, xvar, partial_residuals, bounds, normalizer, multiplier, type = "original", ynew) {

  prob <- integrate(
    density_function, lower = bounds[1], upper = beta,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    normalizer = normalizer, multiplier = multiplier, type = type, ynew = ynew
  )$value

  return(p - prob)

}

post_quant <- function(sig, post_mode, sigma, rate, xvar, partial_residuals, type = "original", ynew) {

  ## Normalizing constant so that the log posterior = 0 at the posterior mode
  multiplier <- -density_function(
    x = post_mode, rate = rate, partial_residuals = partial_residuals,
    sigma = sigma, xvar = xvar, type = type, log = TRUE, ynew = ynew
  )

  ## Determine the normalizing constant (so now the posterior integrates to 1)
  denom <- integrate(
    density_function, lower = -Inf, upper = Inf,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma,
    xvar = xvar, multiplier = multiplier, type = type, ynew = ynew
  )$value

  ## Find the largest density, used for determining bounds for uniroot
  ymode <- density_function(
    x = post_mode, rate = rate, partial_residuals = partial_residuals,
    sigma = sigma, xvar = xvar,
    normalizer = denom,
    multiplier = multiplier,
    type = type, ynew = ynew
  )

  ## Determine bounds
  step <- sigma
  curr <- step
  while (TRUE) {

    xvals <- post_mode + c(-1, 1)*curr
    yvals <- density_function(
      x = xvals, rate = rate, partial_residuals = partial_residuals,
      sigma = sigma, xvar = xvar,
      normalizer = denom,
      multiplier = multiplier,
      type = type, ynew = ynew
    )

    if (all(yvals < (ymode / 1000))) {
      break
    } else {
      curr <- curr + step
    }
  }

  ## redetermine
  denom <- integrate(
    density_function, lower = post_mode - curr, upper = post_mode + curr,
    rate = rate, partial_residuals = partial_residuals, sigma = sigma, xvar = xvar,
    multiplier = multiplier, type = type, ynew = ynew
  )$value

  ## Determine the bounds
  p_lower <- (1 - sig) / 2
  p_upper <- sig + p_lower

  lower <- uniroot(
    obj, c(post_mode - curr, post_mode + curr), p = p_lower, sigma = sigma,
    rate = rate, xvar = xvar, partial_residuals = partial_residuals,
    bounds = c(post_mode - curr, post_mode + curr), normalizer = denom,
    multiplier = multiplier, type = type, ynew = ynew
  )
  upper <- uniroot(
    obj, c(post_mode - curr, post_mode + curr), p = p_upper, sigma = sigma,
    rate = rate, xvar = xvar, partial_residuals = partial_residuals,
    bounds = c(post_mode - curr, post_mode + curr), normalizer = denom,
    multiplier = multiplier, type = type, ynew = ynew
  )
  return(c(lower$root, upper$root))

}

plot_boot <- function(eb_boot, n = 30) {

  rm_lower <- apply(eb_boot[["lower"]], 2, function(x) sum(is.na(x))); names(rm_lower) <- names(eb_boot$truth)
  rm_upper <- apply(eb_boot[["upper"]], 2, function(x) sum(is.na(x))); names(rm_upper) <- names(eb_boot$truth)

  rm_lower <- rm_lower[rm_lower != 0]
  rm_upper <- rm_upper[rm_upper != 0]

  if (length(rm_lower) > 0 ) {
    print(paste0(length(rm_lower), " total variables with NA entries, summary: "))
    print(summary(rm_lower))
  }

  lowers <- apply(eb_boot[["lower"]], 2, mean, na.rm = TRUE)
  uppers <- apply(eb_boot[["upper"]], 2, mean, na.rm = TRUE)
  plot_res <- data.frame(truth = eb_boot[["truth"]], grp = names(eb_boot[["truth"]]), lower = lowers, upper = uppers) %>%
    dplyr::arrange(desc(abs(truth))) %>%
    head(n)

  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = grp)) +
    geom_point(aes(x = truth, y = grp)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")
}

dens <- function(x, z, lambda, sigma2, n, normalizer = 1, multiplier = 1) {

  return(((exp(-(n/sigma2)*(.5*(x - z)^2 + lambda*abs(x)))) / multiplier) / normalizer)

}

obj_simp <- function(beta, p, z, lambda, sigma2, n, normalizer, multiplier, lwr) {

  prob <- integrate(
    dens, lower = lwr, upper = beta,
    z = z, lambda = lambda, sigma2 = sigma2, n = n,
    multiplier = multiplier, normalizer = normalizer
  )$value

  return(p - prob)

}

eb_boot <- function(beta, p = 60, b = 2, n = 100, nboot = 100, type = "original", prog = FALSE, sgm = 1, debias = FALSE, dat = NULL) {

  nlambda <- 100

  if (is.null(dat)) {

    dat <- genDataABN(beta = beta, p = p, a = length(beta), b = b, n = n, sgm = sgm)

    X <- dat$X
    y <- dat$y
    p <- ncol(X)

    cv_res <- cv.ncvreg(X, y, penalty = "lasso")
    sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
    lam <- cv_res$lambda.min
    rate <- (lam*n / sigma2)

    tbeta <- dat$beta

  } else {

    X <- dat$X
    y <- dat$y
    p <- ncol(X)

    cv_res <- cv.ncvreg(X, y, penalty = "lasso")
    sigma2 <- cv_res$cve[cv_res$lambda == cv_res$lambda.min]; sigma <- sqrt(sigma2)
    lam <- cv_res$lambda.min
    rate <- (lam*n / sigma2)

    tbeta <- coef(cv_res$fit, lambda = lam)[-1]

  }

  lowers <- matrix(nrow = nboot, ncol = p)
  tmp_lower <- matrix(nrow = nboot, ncol = p)
  uppers <- matrix(nrow = nboot, ncol = p)

  if (prog) pb <- txtProgressBar(1, nboot, style=3)

  for (i in 1:nboot) {

    idx_new <- sample(1:length(y), replace = TRUE)
    ynew <- y[idx_new]
    xnew <- X[idx_new,,drop=FALSE]
    xnew <- ncvreg::std(xnew)

    lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
    lambda_min <- lam - lam / 100 ## set min to be slightly smaller
    if (lambda_min > lambda_max | lam > lambda_max) {
      lambda_max <- lam + lam / 100
      nlambda <- 2
    }
    lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))

    lasso_fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
    coefs <- coef(lasso_fit, lambda = lam)

    ## Beta specific
    for (j in attr(xnew, "nonsingular")) {

      idx <- as.numeric(which(attr(xnew, "nonsingular") == j))
      post_mode <- coefs[-1][idx]
      xvar <- xnew[,idx,drop=FALSE]

      partial_residuals <- ynew - (coefs[1] + xnew[,-idx,drop=FALSE] %*% coefs[-1][-idx])

      ## Need to compress this into a single call in future
      if (type == "univariate") {

        z <- (1/length(partial_residuals))*as.numeric(t(xvar) %*% partial_residuals)
        R <-  (1/length(partial_residuals))*as.numeric(t(partial_residuals) %*% partial_residuals)

        lwr2 <- pnorm(0, z + lam, sqrt(sigma2 / n))*sqrt(2*pi*sigma2)^(-n) * ((n*lam) / (2*sigma2)) * exp(-(n/(2*sigma2))*((R) - (z + lam)^2))*(sqrt(2*pi*(sigma2/n)))
        upr2 <- pnorm(0, z - lam, sqrt(sigma2 / n), lower.tail = FALSE)*sqrt(2*pi*sigma2)^(-n) * ((n*lam) / (2*sigma2)) * exp(-(n/(2*sigma2))*((R) - (z - lam)^2))*(sqrt(2*pi*(sigma2/n)))

        tdens <- lwr2 + upr2
        prop_lw <- lwr2  / tdens
        prop_up <- upr2  / tdens
        obs_lw <- pnorm(0, z + lam, sqrt(sigma2 / n))
        obs_up <- pnorm(0, z - lam, sqrt(sigma2 / n), lower.tail = FALSE)

        max_lower <- ((pnorm(0, z + lam, sqrt(sigma2 / n)) / obs_lw) * prop_lw)

        tmp <- qnorm(.1 * (obs_lw / prop_lw), z + lam, sqrt(sigma2 / n))

        multiplier <- dens(
          x = post_mode, z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
        )

        denom <- integrate(
          dens, lower = -Inf, upper = Inf,
          z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
          multiplier = multiplier
        )$value


        ymode <- dens(
          x = post_mode, z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
          normalizer = denom, multiplier = multiplier
        )

        step <- 1
        curr <- step
        while (TRUE) {

          xvals <- post_mode + c(-1, 1)*curr
          yvals <- dens(
            x = xvals, z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
            normalizer = denom, multiplier = multiplier
          )

          if (all(yvals < (ymode / 100))) {
            break
          } else {
            curr <- curr + step
          }
        }

        # denom <- integrate(
        #   dens, lower = post_mode - curr, upper = post_mode + curr,
        #   z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
        #   multiplier = multiplier
        # )$value
        #
        # lower <- uniroot(
        #   obj_simp, c(post_mode - curr, post_mode + curr), p = .1,
        #   z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
        #   normalizer = denom, multiplier = multiplier, lwr = post_mode - curr
        # )$root
        # upper <- uniroot(
        #   obj_simp, c(post_mode - curr, post_mode + curr), p = .9,
        #   z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
        #   normalizer = denom, multiplier = multiplier, lwr = post_mode - curr
        # )$root

        lower <- uniroot(
          obj_simp, c(post_mode - curr, post_mode + curr), p = .1,
          z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
          normalizer = denom, multiplier = multiplier, lwr = -Inf
        )$root
        upper <- uniroot(
          obj_simp, c(post_mode - curr, post_mode + curr), p = .9,
          z = z, lambda = lam, sigma2 = sigma2, n = length(partial_residuals),
          normalizer = denom, multiplier = multiplier, lwr = -Inf
        )$root

        # print((prop_lw*obs))
        # print(.1 / pnorm(lower, z - lam, sqrt(sigma2 / n)))

        bounds <- (c(lower, upper) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[j])^(-1)
        tmp_bound_lower <- tmp * (attr(xnew, "scale")[j])^(-1)
        if (tmp_bound_lower > 0) print("Fix me!!!")

      } else if (type == "original") {

        bounds <- (post_quant(.8, post_mode, sigma, rate, xvar, partial_residuals, type = type, ynew) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[idx])^(-1)

      } else if (type == "normal") {

        n <- length(ynew)
        norm_mean <- (2*sigma2*t(xvar) %*% partial_residuals) * (2*sigma2*t(xvar) %*% xvar + n^2*lam^2)^(-1)
        norm_var <- (2*sigma2^2) / (2*(t(xvar) %*% xvar)*sigma2 + n^2*lam^2)
        bounds <- (qnorm(c(.1, .9), norm_mean, sqrt(norm_var)) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[j])^(-1)

      } else if (type == "cadillac") {

        t2i_scale <- lam^2
        r <- ynew - (coefs[1] + xnew[,-idx] %*% coefs[-1][-idx])
        beta <- coefs[-1][idx]
        tau2i_mu <- sqrt((sigma2*lam^2) / beta^2)
        tau2 <- 1 / rinvgauss(1, tau2i_mu, t2i_scale)

        A <- (t(xnew[,idx, drop=FALSE]) %*% xnew[,idx,drop=FALSE]) + (1/tau2)
        mu <- solve(A)*(t(xnew[,idx,drop=FALSE]) %*% r)
        bounds <- (qnorm(c(.1, .9), mu, sqrt(sigma2*solve(A))) + sign(post_mode)*debias*lam*(abs(post_mode) > lam)) * (attr(xnew, "scale")[j])^(-1)

      } else {

        stop(paste0("Type: ", type, " not an option."))

      }

      uppers[i,j] <- bounds[2]
      lowers[i,j] <- bounds[1]
      tmp_lower[i,j] <- tmp_bound_lower

    }

    if(prog) setTxtProgressBar(pb, i)

  }

  print(mean(tmp_lower))
  return(list("lower" = lowers, "upper" = uppers, "truth" = tbeta))

}

eb_boot_sim <- function(beta, p = 60, b = 2, n = 100, nboot = 100, nsim = 100, type = "original", debias = FALSE) {

  overall_cov <- numeric(nsim)
  indiv_cov <- matrix(nrow = nsim, ncol = p)

  pb <- txtProgressBar(1, nsim, style=3)

  ## p, beta, n
  for (iter in 1:nsim) {

    res <- eb_boot(beta = beta, p = p, b = b, n = n, nboot = nboot, type = type, debias = debias)

    lower_int <- apply(res[["lower"]], 2, mean)
    upper_int <- apply(res[["upper"]], 2, mean)
    tbeta <- res[["truth"]]

    indiv_cov[iter,] <- tbeta >= lower_int & tbeta <= upper_int
    overall_cov[iter] <- mean(indiv_cov[iter,])

    setTxtProgressBar(pb, iter)

  }

  return(list("overall_cov" = overall_cov, "indiv_cov" = indiv_cov, "truth" = tbeta))

}

create_plots <- function(tmp) {

  plot_data <- tmp$plots

  # Initialize a list to store the plots
  plotly_plots <- list()

  # Iterate through the plot_data list
  for (i in seq_along(plot_data)) {
    # Extract the data for this plot
    data <- plot_data[[i]][[1]]

    # Create a data frame for ggplot
    plot_df <- data.frame(
      x = data$x,
      posterior = data$posterior,
      prior = data$prior,
      lik = data$lik
    )

    # Create the ggplot object
    p <- ggplot(plot_df, aes(x = x)) +
      geom_line(aes(y = posterior), color = "purple") +
      geom_line(aes(y = prior), color = "red") +
      geom_line(aes(y = lik), color = "blue") +
      geom_vline(xintercept = plot_data[[i]]$post_mode, color = "purple") +
      geom_vline(xintercept = plot_data[[i]]$lower, color = "purple", linetype="dashed") +
      geom_vline(xintercept = plot_data[[i]]$upper, color = "purple", linetype="dashed") +
      labs(title = paste("Lambda: ", round(plot_data[[i]]$lambda, 4)), x = "Beta", y = "Posterior Density")

    # Convert to ggplotly and store in the list
    plotly_plots[[i]] <- ggplotly(p)

  }

  return(plotly_plots)

}
