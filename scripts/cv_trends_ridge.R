########################
## Setup ###############
########################
rm(list = ls())

library(progress) ## install.packages("progress")
library(glue)
library(ncvreg)
library(dplyr)
library(ggplot2)

iterations <- 1000
ridge_loocv_sses <- ridge_2_cv_sses <- ridge_5_cv_sses <- ridge_10_cv_sses <- numeric(iterations)
ridge_loocv_mses <- ridge_2_cv_mses <- ridge_5_cv_mses <- ridge_10_cv_mses <- numeric(iterations)

## Data Scenarios
data_types <- list(
  list(p = 20, beta = c(2, -2, rep(0, 18))),
  list(p = 100, beta = c(1, -1, rep(0, 98))),
  list(p = 50, beta = c(rep(0.5, 8) * rep(c(-1, 1), each = 4), rep(0, 42))),
  list(p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60)))
  # list(p = 100, beta = c(rep(0.25, 32) * rep(c(-1, 1), each = 16), rep(0, 68)))
)
data_type <- 1

res <- list()
for (k in 1:length(data_types)) {
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = iterations, clear = FALSE, width= 100
  )
  for (i in 1:iterations) {

    ##########
    ## Data ##
    ##########
    data <- hdrm::gen_data(n = 50, p = data_types[[data_type]]$p, beta = data_types[[data_type]]$beta)

    ridge_loocv <- cv_ridge(data$X, data$y, nfolds = 50)
    ridge_loocv_beta <- ridge_loocv$fit$beta[-1,ridge_loocv$min]
    ridge_loocv_mse <- sum((data$beta - ridge_loocv_beta)^2)

    ##############
    ## Using CV ##
    ##############
    ridge_10_cv <- cv_ridge(data$X, data$y, nfolds = 10)
    ridge_10_cv_beta <- ridge_10_cv$fit$beta[-1,ridge_10_cv$min]
    ridge_5_cv <- cv_ridge(data$X, data$y, nfolds = 5)
    ridge_5_cv_beta <- ridge_5_cv$fit$beta[-1,ridge_5_cv$min]
    ridge_2_cv <- cv_ridge(data$X, data$y, nfolds = 2)
    ridge_2_cv_beta <- ridge_2_cv$fit$beta[-1,ridge_2_cv$min]

    ####################
    ## Calculate SSEs ##
    ####################
    ridge_10_cv_mses[i] <- sum((data$beta - ridge_10_cv_beta)^2) - ridge_loocv_mse
    ridge_5_cv_mses[i] <- sum((data$beta - ridge_5_cv_beta)^2) - ridge_loocv_mse
    ridge_2_cv_mses[i] <- sum((data$beta - ridge_2_cv_beta)^2) - ridge_loocv_mse

    pb$tick()

  }
  res[[k]] <- tibble(
    folds = c(rep(10, iterations), rep(5, iterations), rep(2, iterations)),
    sse = c(ridge_10_cv_mses, ridge_5_cv_mses, ridge_2_cv_mses)
  ) %>%
    mutate(scenario = glue("Scenario {k}"))

}

indexr::save_objects("./rds", bind_rows(res), args_list = list(script_name = "cv_trends_ridge"), overwrite = TRUE)

