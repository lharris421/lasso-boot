########################
## Setup ###############
########################
rm(list = ls())

library(progress) ## install.packages("progress")
library(glue)
library(ncvreg)
library(hdrm)

iterations <- 100
lambdas_1se <- lambdas_min <- lambdas_oracle <- numeric(iterations)

## Data Scenarios
data_types <- list(
  list(n = 50, p = 20, beta = c(2, -2, rep(0, 18))),
  list(n = 50, p = 100, beta = c(1, -1, rep(0, 98))),
  list(n = 50, p = 50, beta = c(rep(0.5, 8) * rep(c(-1, 1), each = 4), rep(0, 42))),
  list(n = 50, p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60))),
  list(n = 100, p = 100, beta = c(2, -2, rep(0, 98)), corr = "autoregressive", rho = 0.7)
)
data_type <- 5

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
  cv_fit <- cv.ncvreg(data$X, data$y, penalty = "lasso")
  lambdas_min[i] <- cv_fit$lambda.min
  lambdas_1se[i] <- cv_fit$lambda[min(which(cv_fit$cve <= cv_fit$cve[cv_fit$min] + cv_fit$cvse[cv_fit$min]))]

  ###################
  ## Oracle lambda ##
  ###################
  lambdas_oracle[i] <- cv_fit$lambda[which.min(colMeans((data$beta - cv_fit$fit$beta[-1,])^2))]

  pb$tick()

}

mean(lambdas_min)
mean(lambdas_oracle)
mean(lambdas_1se)

hist(lambdas_min - lambdas_oracle)

sd(lambdas_min)
sd(lambdas_oracle)
