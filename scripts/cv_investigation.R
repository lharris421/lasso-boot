########################
## Setup ###############
########################
rm(list = ls())
source("./scripts/setup.R")

library(progress) ## install.packages("progress")
library(glue)
library(ncvreg)
library(hdrm)

iterations <- 100
lambdas_pic <- lambdas_new <- lambdas_aic <- lambdas_bic <- lambdas_min <- lambdas_oracle <- lambdas_oracle2 <- numeric(iterations)

## Data Scenarios
data_types <- list(
  list(n = 50, p = 20, beta = c(2, -2, rep(0, 18))),
  list(n = 50, p = 100, beta = c(1, -1, rep(0, 98))),
  list(n = 50, p = 50, beta = c(rep(0.5, 8) * rep(c(-1, 1), each = 4), rep(0, 42))),
  list(n = 50, p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60))),
  list(n = 50, p = 100, beta = c(2, -2, rep(0, 98)), corr = "autoregressive", rho = 0.7)
)
data_type <- 2

pb <- progress_bar$new(
  format = "[:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 100
)
for (i in 1:iterations) {

  ##########
  ## Data ##
  ##########
  data <- do.call("gen_data", data_types[[data_type]])

  ##############
  ## Using CV ##
  ##############
  cv_fit <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 50)
  lambdas_min[i] <- cv_fit$lambda.min
  lambdas_aic[i] <- cv_fit$lambda[which.min(AIC(cv_fit$fit))]
  lambdas_bic[i] <- cv_fit$lambda[which.min(BIC(cv_fit$fit))]
  lambdas_pic[i] <- cv_fit$lambda[which.min(pic(data$X, data$y, relaxed = TRUE))]


  ###################
  ## Oracle lambda ##
  ###################
  # lambdas_oracle2[i] <- cv_fit$lambda[which.min(colMeans((data$y - cbind(1, data$X) %*% cv_fit$fit$beta)^2))]
  lambdas_oracle[i] <- cv_fit$lambda[which.min(colMeans((data$beta - cv_fit$fit$beta[-1,])^2))]

  pb$tick()

}

mean(lambdas_oracle)
mean(lambdas_min)
mean(lambdas_aic)
mean(lambdas_bic)
mean(lambdas_pic)

