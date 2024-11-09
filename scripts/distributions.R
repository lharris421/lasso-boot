## Setup
source("./scripts/setup.R")

params <- list(seed = 1234, iterations = 1,
               simulation_function = "gen_data_distribution", simulation_arguments = list(
               ), script_name = "distributions")


for (i in 1:length(methods)) {
  methods[[i]]$method_arguments["alpha"] <- 0.2
}
methods <- methods[c("lasso_proj_boot")]

ns <- NULL
rhos <- NULL
correlations <- NULL
distributions <- c("Scheetz2006")
true_lambda <- NULL
true_sigma2 <- NULL

simulations <- expand.grid(
  # "n" = ns,
  # "rho" = rhos,
  # "corr" = correlations,
  "distribution" = distributions,
  #"true_lambda" = true_lambda,
  #"true_sigma2" = true_sigma2,
  stringsAsFactors = FALSE
)

pb <- txtProgressBar(min = 0, max = params$iterations, initial = 0, style = 3)

set.seed(params$seed)
seeds <- round(runif(params$iterations) * 1e9)
for(k in 1:nrow(simulations)) {

  ## Need to make more flexible
  if (!is.null(ns)) params$simulation_arguments$n <- simulations %>% filter(row_number() == k) %>% pull(n)
  if (!is.null(distributions)) params$simulation_arguments$distribution <- simulations %>% filter(row_number() == k) %>% pull(distribution)
  if (!is.null(rhos)) params$simulation_arguments$rho <- simulations %>% filter(row_number() == k) %>% pull(rho)
  if (!is.null(correlations)) params$simulation_arguments$corr <- simulations %>% filter(row_number() == k) %>% pull(corr)
  if (!is.null(true_lambda)) params$true_lambda <- simulations %>% filter(row_number() == k) %>% pull(true_lambda)
  if (!is.null(true_sigma2)) params$true_sigma2 <- simulations %>% filter(row_number() == k) %>% pull(true_sigma2)

  print(params)

  res <- list()
  for(i in 1:params$iterations) {

    set.seed(seeds[i])
    data <- do.call(params$simulation_function, params$simulation_arguments)
    lasso_lambda <- NULL
    lasso_sigma2 <- NULL
    mcp_lambda <- NULL
    mcp_sigma2 <- NULL
    if (!is.null(params$same_lambda) && params$same_lambda) {
      glmnet_cv <- cv.glmnet(data$X, data$y)
      lasso_lambda <- glmnet_cv$lambda.min
    } else if (!is.null(params$true_lambda) && params$true_lambda) {
      lasso_lambda <- sqrt(2*params$simulation_arguments$p) / params$simulation_arguments$n
    } else if (any(c("lasso", "lasso_relaxed", "lasso_boot") %in% names(methods))) {
      lasso_cv <- cv.ncvreg(data$X, data$y, penalty = "lasso")
      lasso_lambda <- lasso_cv$lambda.min
      lasso_sigma2 <- lasso_cv$cve[lasso_cv$min]
    }
    if (any(c("mcp", "mcp_relaxed", "mcp_boot") %in% names(methods))) {
      mcp_cv <- cv.ncvreg(data$X, data$y, penalty = "MCP")
      mcp_lambda <- mcp_cv$lambda.min
      mcp_sigma2 <- mcp_cv$cve[mcp_cv$min]
    }
    if (!is.null(params$true_sigma2) && params$true_sigma2) {
      lasso_sigma2 <- 1
    }

    res[[i]] <- list()
    for (j in 1:length(methods)) {

      if (!is.null(params$same_lambda) && params$same_lambda) {
        lambda <- lasso_lambda
      } else if (names(methods)[j] %in% c("mcp", "mcp_relaxed", "mcp_boot")) {
        lambda <- mcp_lambda
      } else if (names(methods)[j] %in% c("lasso", "lasso_relaxed", "lasso_boot")) {
        lambda <- lasso_lambda
      } else {
        lambda <- NULL
      }

      if (names(methods[j]) == "lasso_boot" & !is.null(lasso_sigma2)) {
        time_taken <- system.time({
          results <- do.call(methods[[j]]$method, c(list(X = data$X, y = data$y, lambda = lambda, sigma2 = lasso_sigma2), methods[[j]]$method_arguments))
        })
      } else if (names(methods[j]) == "mcp_boot") {
        time_taken <- system.time({
          results <- do.call(methods[[j]]$method, c(list(X = data$X, y = data$y, lambda = lambda, sigma2 = mcp_sigma2), methods[[j]]$method_arguments))
        })
      } else if (!is.null(lambda)) {
        time_taken <- system.time({
          results <- do.call(methods[[j]]$method, c(list(X = data$X, y = data$y, lambda = lambda), methods[[j]]$method_arguments))
        })
      } else {
        time_taken <- system.time({
          results <- do.call(methods[[j]]$method, c(list(X = data$X, y = data$y), methods[[j]]$method_arguments))
        })
      }

      res[[i]][[names(methods)[j]]] <- results %>%
        mutate(truth = data$beta, iteration = i, method = names(methods)[j], time = time_taken["elapsed"])

    }

    setTxtProgressBar(pb,i)

  }

  method_names <- names(methods)
  for (j in 1:length(methods)) {

    ## Combine results for a give method
    results <- bind_rows(extract_named_elements(res, method_names[j]))

    ## Add method info to params
    method_info <- methods[[method_names[j]]]
    save_params <- params
    save_params[names(method_info)] <- method_info

    ## Save results
    indexr::save_objects("./rds", results, args_list = save_params, overwrite = TRUE)

  }

}
