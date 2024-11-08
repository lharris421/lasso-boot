source("./scripts/setup.R")

## A script version set up of the simulation
## this differs mildly in behavior, allowing a more granular handling in lambda

params <- list(seed = 1234, iterations = 1000,
               simulation_function = "gen_data_abn", simulation_arguments = list(
                 n = 100, p = 100, a = 1, b = 1, rho = 0.5,
                 beta = 2 ## Should jsut enter the non zeros
               ), script_name = "b_selection")

res <- list()

pb <- txtProgressBar(min = 0, max = params$iterations, initial = 0, style = 2) 

set.seed(params$seed)
seeds <- round(runif(params$iterations) * 1e9)
origDebiasedEstAs <- origABCorrs <- origANBiases <- origErrBiases <- orig_ns <- orig_bs <- orig_as <- numeric(1000)
Ns <- As <- Bs <- debiasedEstAs <- abCorrs <- anBiases <- ErrBiases <- matrix(nrow = 1000, ncol = 1000)
for(i in 1:params$iterations) {
  
  set.seed(seeds[i])
  data <- do.call(params$simulation_function, params$simulation_arguments)
  cv_fit <- cv.ncvreg(data$X, data$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min / max(cv_fit$fit$lambda)
  errs <- data$y - (data$X %*% data$beta)
  if (abs(lambda - 0.05) < 1e-9) {lambda <- 0.050001}
  if (abs(lambda - 1) < 1e-9) {lambda <- 0.99999}
  
  orig_coefs <- coef(cv_fit$fit, lambda = max(cv_fit$fit$lambda) * lambda)[-1]
  orig_as[i] <- orig_coefs[1]
  orig_bs[i] <- orig_coefs[2]
  orig_ns[i] <- sum(orig_coefs[3:100])
  origStdX <- ncvreg::std(data$X)
  origABCorrs[i] <- (1/100)* origStdX[,1] %*% origStdX[,2]
  origANBiases[i] <- sum(((1/100)* origStdX[,1] %*% origStdX[,3:100]) * -orig_coefs[3:100])
  origErrBiases[i] <- (1/100) * origStdX[,1] %*% errs * (attr(origStdX, "scale")^(-1))[1]
  
  y <- data$y - mean(data$y)
  origStdFit <- ncvreg(origStdX, y, penalty = "lasso")
  origModes <- coef(origStdFit, lambda = max(origStdFit$lambda) * lambda)[-1]
  orig_partial_residuals <- y - (
    as.numeric(origStdX %*% origModes) - (origStdX * matrix(origModes, nrow = nrow(origStdX), ncol = ncol(origStdX), byrow=TRUE))
  )
  orig_z <- (1/100)*colSums(origStdX * orig_partial_residuals)
  orig_z <- orig_z * attr(origStdX, "scale")^(-1)
  origDebiasedEstAs[i] <- orig_z[1]

  n <- a <- b <- debiasedEstA <- abCorr <- anBias <- ErrBias <- numeric(1000)
  anAddCorrs <- matrix(nrow = 1000, ncol = 98)
  anAddBias <- matrix(nrow = 1000, ncol = 98)
  for (j in 1:1000) {
    
    boot_sample <- sample(1:100, replace = TRUE)
    xnew <- data$X[boot_sample,]
    ynew <- data$y[boot_sample]
    ynew <- ynew - mean(ynew)
    stdX <- std(xnew)
    
    fit <- ncvreg(xnew, ynew, penalty = "lasso")
    stdFit <- ncvreg(stdX, ynew, penalty = "lasso")
    coefs <- coef(fit, lambda = lambda * max(fit$lambda))[-1]
    
    n[j] <- sum(coefs[3:100])
    b[j] <- coefs[2]
    a[j] <- coefs[1]
    
    abCorr[j] <- (1/100) * stdX[,1] %*% stdX[,2]
    anCorrs <- (1/100) * stdX[,1] %*% stdX[,3:100]
    anBias[j] <- sum(anCorrs * -coefs[3:100])
    ErrBias[j] <- (1/100) * stdX[,1] %*% errs[boot_sample] * (attr(stdX, "scale")^(-1))[1]
    
    anAddBias[j,] <- coefs[3:100] - orig_coefs[3:100]
    anAddCorrs[j,] <- anCorrs - ((1/100)* origStdX[,1] %*% origStdX[,3:100])
    
    modes <- coef(stdFit, lambda = max(stdFit$lambda) * lambda)[-1]
    partial_residuals <- ynew - (
      as.numeric(stdX %*% modes) - (stdX * matrix(modes, nrow = nrow(stdX), ncol = ncol(stdX), byrow=TRUE))
    )
    z <- (1/100)*colSums(stdX * partial_residuals)
    z <- z * attr(stdX, "scale")^(-1)
    debiasedEstA[j] <- z[1]
    
  }
  Ns[i,] <- n
  As[i,] <- a
  Bs[i,] <- b
  debiasedEstAs[i,] <- debiasedEstA
  abCorrs[i,] <- abCorr
  anBiases[i,] <- anBias
  ErrBiases[i,] <- ErrBias
  
  setTxtProgressBar(pb,i)
  
}

results <- list(
  "orig_est" = origDebiasedEstAs, "orig_corrs" = origABCorrs,
  "an_bias" = anBiases, "err_bias" = ErrBiases,
  "orig_an_bias" = origANBiases, "orig_err_bias" = origErrBiases,
  "orig_ns" = orig_ns, "orig_bs" = orig_bs, "orig_as" = orig_as,
  "ests" = debiasedEstAs, "corrs" = abCorrs,
  "ns" = Ns, "bs" = Bs, "as" = As
)
indexr::save_objects("./rds", results, args_list = params, overwrite = TRUE)
  