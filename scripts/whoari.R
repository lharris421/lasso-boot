source("./scripts/setup.R")

dat <- hdrm::readData("whoari")

cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
cv_fit_mcp <- cv.ncvreg(dat$X, dat$y, penalty = "MCP")

for (i in 1:length(methods)) {
  curr_lambda <- ifelse(methods[[i]]$method == "mcp", cv_fit_mcp$lambda.min, cv_fit$lambda.min)
  res <- do.call(methods[[i]]$method, c(list(X = dat$X, y = dat$y, lambda = curr_lambda), methods[[i]]$method_arguments))
  indexr::save_objects("./rds", res, args_list = methods[[i]], overwrite = TRUE)
}

