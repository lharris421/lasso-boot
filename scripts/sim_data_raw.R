source("./scripts/setup.R")

params <- list(seed = 1234, iterations = 1000,
               simulation_function = "gen_data", simulation_arguments = list(
                 n = 100, p = 100,
                 beta = c(-2, 2, -1, 1, -0.5, 0.5, -0.5, 0.5, rep(0, 92))
               ), script_name = "sim_data_raw")

res <- list()

set.seed(params$seed)
data <- do.call(params$simulation_function, params$simulation_arguments)
lambda <- cv.ncvreg(data$X, data$y, penalty = "lasso")$lambda.min

res_pipe <- pipe_ncvreg(X = data$X, y = data$y, lambda = lambda)
res_debiased <- boot_ncv(X = data$X, y = data$y, lambda = lambda, return_boot = TRUE)

res <- list("pipe" = res_pipe, "debiased" = res_debiased)

indexr::save_objects("./rds", res, args_list = params, overwrite = TRUE)


