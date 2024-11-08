source("./scripts/setup.R")

params <- list(
  seed = 1234,
  iterations = 1000,
  simulation_function = "gen_data_abn",
  simulation_arguments = list(
    n = 100, p = 100, a = 1, b = 1, rho = 0.99,
    beta = 2
  ),
  script_name = "bias_decomposition"
)

bootstrap_samples <- 1000
res <- list()

pb <- txtProgressBar(min = 0, max = params$iterations, initial = 0, style = 2)

set.seed(params$seed)
seeds <- round(runif(params$iterations) * 1e9)

lambdas <- numeric(params$iterations)
for (i in 1:params$iterations) {
  
  set.seed(seeds[i])
  data <- do.call(params$simulation_function, params$simulation_arguments)
  cv_fit <- cv.ncvreg(data$X, data$y, penalty = "lasso")
  lambdas[i] <- cv_fit$lambda.min
  lambda <- cv_fit$lambda.min / max(cv_fit$fit$lambda)
  errs <- data$y - (data$X %*% data$beta)
  if (abs(lambda - 0.05) < 1e-9) { lambda <- 0.050001 }
  if (abs(lambda - 1) < 1e-9) { lambda <- 0.99999 }
  
  # Get variable names and indices
  var_names <- names(data$beta)
  a1_index <- which(var_names == "A1")  # Focus on variable "A1"
  other_a_indices <- which(grepl("^A", var_names) & var_names != "A1")  # Other A variables
  b_indices <- which(grepl("^B", var_names))
  n_indices <- which(grepl("^N", var_names))
  n_obs <- nrow(data$X)
  
  # Initialize storage variables on the first iteration
  if (i == 1) {
    orig_debiased <- numeric(params$iterations)
    orig_n_biases <- numeric(params$iterations)
    orig_err_biases <- numeric(params$iterations)
    true_betas <- data$beta
    
    n_other_a <- length(other_a_indices)
    n_b <- length(b_indices)
    if (n_other_a > 0) {
      orig_a_biases <- matrix(nrow = params$iterations, ncol = n_other_a)
      colnames(orig_a_biases) <- var_names[other_a_indices]
    } else {
      orig_a_biases <- NULL
    }
    if (n_b > 0) {
      orig_b_biases <- matrix(nrow = params$iterations, ncol = n_b)
      colnames(orig_b_biases) <- var_names[b_indices]
    } else {
      orig_b_biases <- NULL
    }
    
    debiaseds <- matrix(nrow = params$iterations, ncol = bootstrap_samples)
    n_biases <- matrix(nrow = params$iterations, ncol = bootstrap_samples)
    err_biases <- matrix(nrow = params$iterations, ncol = bootstrap_samples)
    
    if (n_other_a > 0) {
      a_biases <- array(dim = c(params$iterations, bootstrap_samples, n_other_a))
      dimnames(a_biases)[[3]] <- var_names[other_a_indices]
    } else {
      a_biases <- NULL
    }
    if (n_b > 0) {
      b_biases <- array(dim = c(params$iterations, bootstrap_samples, n_b))
      dimnames(b_biases)[[3]] <- var_names[b_indices]
    } else {
      b_biases <- NULL
    }
  }
  
  orig_coefs <- coef(cv_fit$fit, lambda = max(cv_fit$fit$lambda) * lambda)[-1]
  origStdX <- ncvreg::std(data$X)
  
  # Compute correlations between A1 and other non-null predictors
  if (n_other_a > 0) {
    orig_a_biases[i, ] <- (orig_coefs[other_a_indices]-data$beta[other_a_indices]) * as.numeric((1 / n_obs) * t(origStdX[, a1_index]) %*% origStdX[, other_a_indices])
  }
  
  if (n_b > 0) {
    orig_b_biases[i, ] <- (orig_coefs[b_indices]-data$beta[b_indices]) * as.numeric((1 / n_obs) * t(origStdX[, a1_index]) %*% origStdX[, b_indices])
  }
  
  # Compute bias due to N variables for A1
  if (length(n_indices) > 0) {
    an_corrs <- as.numeric((1 / n_obs) * t(origStdX[, a1_index]) %*% origStdX[, n_indices])
    orig_n_biases[i] <- sum(an_corrs * orig_coefs[n_indices])
  } else {
    orig_n_biases[i] <- 0
  }
  
  # Compute error bias for A1
  orig_err_biases[i] <- -(1 / n_obs) * sum(origStdX[, a1_index] * errs) * (attr(origStdX, "scale")^(-1))[a1_index]
  
  # Debiased estimate for A1
  y <- data$y - mean(data$y)
  origStdFit <- ncvreg(origStdX, y, penalty = "lasso")
  origModes <- coef(origStdFit, lambda = max(origStdFit$lambda) * lambda)[-1]
  orig_partial_residuals <- y - (
    as.numeric(origStdX %*% origModes) - (origStdX * matrix(origModes, nrow = nrow(origStdX), ncol = ncol(origStdX), byrow = TRUE))
  )
  orig_z <- (1 / n_obs) * colSums(origStdX * orig_partial_residuals)
  orig_z <- orig_z * attr(origStdX, "scale")^(-1)
  orig_debiased[i] <- orig_z[a1_index]
  
  # Initialize storage for bootstrap samples
  debiased <- numeric(bootstrap_samples)
  n_bias <- numeric(bootstrap_samples)
  err_bias <- numeric(bootstrap_samples)
  
  
  if (n_other_a > 0) {
    a_bias <- matrix(nrow = bootstrap_samples, ncol = n_other_a)
  }
  if (n_b > 0) {
    b_bias <- matrix(nrow = bootstrap_samples, ncol = n_b)
  }
  
  for (j in 1:bootstrap_samples) {
    
    boot_sample <- sample(1:n_obs, replace = TRUE)
    xnew <- data$X[boot_sample, ]
    ynew <- data$y[boot_sample]
    ynew <- ynew - mean(ynew)
    stdX <- std(xnew)
    
    fit <- ncvreg(xnew, ynew, penalty = "lasso")
    stdFit <- ncvreg(stdX, ynew, penalty = "lasso")
    coefs <- coef(fit, lambda = lambda * max(fit$lambda))[-1]
    
    if (n_other_a > 0) {
      a_bias[j, ] <- (coefs[other_a_indices] - data$beta[other_a_indices]) * as.numeric((1 / n_obs) * t(stdX[, a1_index]) %*% stdX[, other_a_indices])
    }
    if (n_b > 0) {
      b_bias[j, ] <- (coefs[b_indices]- data$beta[b_indices]) * as.numeric((1 / n_obs) * t(stdX[, a1_index]) %*% stdX[, b_indices])
    }
    
    # Bias due to N variables for A1
    if (length(n_indices) > 0) {
      n_bias[j] <- sum(as.numeric((1 / n_obs) * t(stdX[, a1_index]) %*% stdX[, n_indices]) * coefs[n_indices])
    } else {
      n_bias[j] <- 0
    }
    
    # Error bias for A1
    err_bias[j] <- -(1 / n_obs) * sum(stdX[, a1_index] * errs[boot_sample]) * (attr(stdX, "scale")^(-1))[a1_index]
    
    # Debiased estimate for A1
    modes <- coef(stdFit, lambda = max(stdFit$lambda) * lambda)[-1]
    partial_residuals <- ynew - (
      as.numeric(stdX %*% modes) - (stdX * matrix(modes, nrow = nrow(stdX), ncol = ncol(stdX), byrow = TRUE))
    )
    z <- (1 / n_obs) * colSums(stdX * partial_residuals)
    z <- z * attr(stdX, "scale")^(-1)
    debiased[j] <- z[a1_index]
    
  }
  debiaseds[i, ] <- debiased
  n_biases[i, ] <- n_bias
  err_biases[i, ] <- err_bias
  
  
  if (n_other_a > 0) {
    a_biases[i, , ] <- a_bias
  }
  if (n_b > 0) {
    b_biases[i, , ] <- b_bias
  }
  
  setTxtProgressBar(pb, i)
  
}

results <- list(
  "orig_est" = orig_debiased,
  "orig_n_bias" = orig_n_biases,
  "orig_a_bias" = orig_a_biases,
  "orig_b_bias" = orig_b_biases,
  "orig_err_bias" = orig_err_biases,
  "ests" = debiaseds,
  "n_bias" = n_biases,
  "a_bias" = a_biases,
  "b_bias" = b_biases,
  "err_bias" = err_biases,
  "lambdas" = lambdas,
  "true_betas" = true_betas
)

indexr::save_objects("./rds", results, args_list = params, overwrite = TRUE)
