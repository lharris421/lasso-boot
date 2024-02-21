source("./scripts/setup/setup.R")

## Parameters
n <- 50
p <- 25
ci_method <- "quantile"
methods <- c("ridge", methods)

## Selected example
selected_example <- sample(1:100, 1)

cis <- list()
examples <- list()
current_seed <- my_seed
for (i in 1:100) {

  print(i)
  current_seed <- current_seed + i
  set.seed(current_seed)
  dat <- gen_data_abn(n = n, p = p, a = 1, b = 1, rho = .99)

  for (j in 1:length(methods)) {
    set.seed(current_seed)
    if (i == 1) {
      cis[[j]] <- list()
    }
    if (methods[j] == "ridge") {
      ### Ridge
      ridge_cv <- cv_ridge(dat$X, dat$y)
      cis[[j]][[i]] <- cbind(confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min), method = "ridge") %>%
        data.frame() %>%
        mutate(group = i)
    } else {
      ### Lasso-boot sample
      lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, max.iter = 1e9, method = methods[j], nboot = nboot)
      cis[[j]][[i]] <- ci.boot.ncvreg(lassoboot, ci_method= ci_method) %>%
        dplyr::mutate(method = methods[j], group = i)
    }

    if (i == selected_example) {
      if (methods[j] == "ridge") {
        examples[[j]] <- cbind(confint(ridge_cv$fit, level = 0.8, lambda = ridge_cv$lambda.min), "Estimate" = coef(ridge_cv$fit, lambda = ridge_cv$lambda.min))
      } else {
        examples[[j]] <- lassoboot
      }
    }

  }

}
names(cis) <- methods
names(examples) <- methods

for (i in 1:length(methods)) {
  confidence_interval <- cis[[methods[i]]]
  example <- examples[[methods[i]]]
  save(confidence_interval, example, file = glue("./rds/highcorr_{methods[i]}.rds"))
}

