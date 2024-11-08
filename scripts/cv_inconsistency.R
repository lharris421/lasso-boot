## Setup
source("./scripts/setup.R")

params <- list(seed = 1234, iterations = 100,
               simulation_function = "gen_data_abn", simulation_arguments = list(
                 n = 100, p = 100, a= 6, b = 2, SNR = 1
               ), script_name = "cv_inconsistency")


pb <- txtProgressBar(min = 0, max = params$iterations, initial = 0, style = 3)

res <- data.frame(simulation = 1:params$iterations)

# set.seed(params$seed)
# seeds <- round(runif(params$iterations) * 1e9)
data <- do.call(params$simulation_function, params$simulation_arguments)

for(i in 1:params$iterations) {
    
  # set.seed(seeds[i]) 
  lasso_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 10)
  res[i,"lambda_max"] <- max(lasso_cv$lambda)
  res[i,"lambda_min"] <- lasso_cv$lambda.min
  res[i,"cve_min"] <- lasso_cv$cve[lasso_cv$min]
  
  # max_allowable_cve <-  lasso_cv$cve[lasso_cv$min] + lasso_cv$cvse[lasso_cv$min]
  # idx_1se <- min(which(lasso_cv$cve <= max_allowable_cve))
  # res[i,"lambda_1se"] <- lasso_cv$lambda[idx_1se]
  # res[i,"cve_1se"] <- lasso_cv$cve[idx_1se]
  
  ncoefs <- apply(coef(lasso_cv$fit)[-1,], 2, function(x) sum(x != 0))
  
  idx_aware <- which.min(lasso_cv$cve + colSums((data$y - lasso_cv$fit$linear.predictors)^2) / (params$simulation_arguments$n - ncoefs))
  res[i,"lambda_aware"] <- lasso_cv$lambda[idx_aware]
  res[i,"cve_aware"] <- lasso_cv$cve[idx_aware]
  
  res[i,"n_coef_min"] <- ncoefs[lasso_cv$min]
  res[i,"n_coef_aware"] <- ncoefs[idx_aware]
  
  res[i,"prediction_error"] <- drop(crossprod(data$y - lasso_cv$fit$linear.predictors[,lasso_cv$min]))
  res[i,"rhee_var"] <- res[i,"prediction_error"] / (params$simulation_arguments$n - res[i,"n_coef_min"])
  

  setTxtProgressBar(pb,i)
  
}
  
## indexr::save_objects("./rds", results, args_list = params, overwrite = TRUE)

table(res$n_coef_min)
table(res$n_coef_aware)

