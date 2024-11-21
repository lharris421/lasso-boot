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
lasso_loocv_sses <- lasso_2_cv_sses <- lasso_5_cv_sses <- lasso_10_cv_sses <- numeric(iterations)
lasso_loocv_mses <- lasso_2_cv_mses <- lasso_5_cv_mses <- lasso_10_cv_mses <- numeric(iterations)

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

    lasso_loocv <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 50)
    lasso_loocv_beta <- coef(lasso_loocv)[-1]
    lasso_loocv_mse <- sum((data$beta - lasso_loocv_beta)^2)

    ##############
    ## Using CV ##
    ##############
    lasso_10_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso")
    lasso_10_cv_beta <- coef(lasso_10_cv)[-1]
    lasso_5_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 5)
    lasso_5_cv_beta <- coef(lasso_5_cv)[-1]
    lasso_2_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 2)
    lasso_2_cv_beta <- coef(lasso_2_cv)[-1]

    ####################
    ## Calculate SSEs ##
    ####################
    lasso_10_cv_mses[i] <- sum((data$beta - lasso_10_cv_beta)^2) - lasso_loocv_mse
    lasso_5_cv_mses[i] <- sum((data$beta - lasso_5_cv_beta)^2) - lasso_loocv_mse
    lasso_2_cv_mses[i] <- sum((data$beta - lasso_2_cv_beta)^2) - lasso_loocv_mse

    pb$tick()

  }
  res[[k]] <- tibble(
    folds = c(rep(10, iterations), rep(5, iterations), rep(2, iterations)),
    sse = c(lasso_10_cv_mses, lasso_5_cv_mses, lasso_2_cv_mses)
  ) %>%
    mutate(scenario = glue("Scenario {k}"))

}

indexr::save_objects("./rds", bind_rows(res), args_list = list(script_name = "cv_trends_lasso"), overwrite = TRUE)

