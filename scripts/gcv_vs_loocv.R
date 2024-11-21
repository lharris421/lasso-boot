## 1) Show that gcv is inferior to loocv
#### 1 hist per simulation with the differences of the sses to show distribution is shiftered from zero
#### and that loocv is in fact better
#### for 1-4

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
  list(p = 100, beta = c(rep(0.3, 40) * rep(c(-1, 1), each = 20), rep(0, 60)))
)

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
    data <- hdrm::gen_data(n = 50, p = data_types[[k]]$p, beta = data_types[[k]]$beta)

    ##############
    ## Using CV ##
    ##############
    ridge_loocv <- cv_ridge(data$X, data$y, nfolds = 50)
    ridge_loocv_beta <- ridge_loocv$fit$beta[-1,ridge_loocv$min]
    ridge_gcv <- hdrm::ridge(data$X, data$y)
    ridge_gcv_beta <-ridge_gcv$beta[-1,which.min(ridge_gcv$GCV)]

    ####################
    ## Calculate SSEs ##
    ####################
    # if (sum((data$beta - ridge_gcv_beta)^2) > 10) stop("gcv is wack")
    ridge_loocv_sses[i] <- sum((data$beta - ridge_loocv_beta)^2)
    ridge_gcv_sses[i] <- sum((data$beta - ridge_gcv_beta)^2)
    pb$tick()

  }
  res[[k]] <- tibble(
    "ridge_gcv" = ridge_gcv_sses,
    "ridge_loocv" = ridge_loocv_sses,
  ) %>%
    mutate(scenario = glue("Scenario {k}"))
}

## Goal 1 %>%
indexr::save_objects("./rds", bind_rows(res), args_list = list(script_name = "gcv_vs_loocv"), overwrite = TRUE)


  mutate(differences = ridge_gcv - ridge_loocv) %>%
  ggplot(aes(x = differences)) +
  geom_histogram(bins = 10) +
  facet_wrap(~scenario, nrow = 4, scales = "free") +
  theme_bw() +
  geom_vline(xintercept = 0, color = "red")
