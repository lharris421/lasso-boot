## 2) Show for simulation 5 how this contributes to lasso looking as good as ridge
#### 3 histograms:
#### a) ridge vs lasso using gcv vs loocv -> shows they are equiv
#### b) ridge vs lasso using loocv vs loocv -> shows not equiv
#### c) ridge vs ridge using gcv vs loocv -> shows not equiv, explains pattern above

## Still need to then bridge gap between loocv and 10 fold cv
## This is goal of other simulation

########################
## Setup ###############
########################
source("./scripts/setup.R")
rm(list = ls())

library(progress) ## install.packages("progress")
library(glue)
library(ncvreg)

iterations <- 1000
lasso_loocv_sses <- ridge_loocv_sses <- ridge_gcv_sses <- numeric(iterations)

## Data Scenarios
data_types <- list(
  list(p = 20, beta = c(2, -2, rep(0, 18))), # p = 25 in HW, betas mag = 1
  list(p = 100, beta = c(1, -1, rep(0, 98))),
  list(p = 50, beta = c(rep(0.5, 8) * rep(c(-1, 1), each = 4), rep(0, 42))),
  list(p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60))),
  list(p = 100, beta = c(rep(0.25, 32) * rep(c(-1, 1), each = 16), rep(0, 68)))
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
  data <- hdrm::gen_data(n = 50, p = data_types[[data_type]]$p, beta = data_types[[data_type]]$beta)

  ##############
  ## Using CV ##
  ##############
  lasso_loocv <- cv.ncvreg(data$X, data$y, penalty = "lasso", nfolds = 50)
  lasso_loocv_beta <- coef(lasso_loocv)[-1]
  ridge_loocv <- cv_ridge(data$X, data$y, nfolds = 50)
  ridge_loocv_beta <- ridge_loocv$fit$beta[-1,ridge_loocv$min]
  ridge_gcv <- hdrm::ridge(data$X, data$y)
  ridge_gcv_beta <-ridge_gcv$beta[-1,which.min(ridge_gcv$GCV)]

  ####################
  ## Calculate SSEs ##
  ####################
  lasso_loocv_sses[i] <- sum((data$beta - lasso_loocv_beta)^2)
  ridge_loocv_sses[i] <- sum((data$beta - ridge_loocv_beta)^2)
  ridge_gcv_sses[i] <- sum((data$beta - ridge_gcv_beta)^2)
  pb$tick()

}

## Goal 2
res <- tibble(
  "ridge_gcv" = ridge_gcv_sses,
  "ridge_loocv" = ridge_loocv_sses,
  "lasso_loocv" = lasso_loocv_sses,
) %>%
  mutate(
    diff_rl_gl = ridge_gcv - lasso_loocv,
    diff_rl_ll = ridge_loocv - lasso_loocv,
    diff_rr_gl = ridge_gcv - ridge_loocv,
  )

indexr::save_objects("./rds", res, args_list = list(script_name = "gcv_incomparability"), overwrite = TRUE)
