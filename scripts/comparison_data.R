source("./scripts/setup/setup.R")
library(tictoc)
# ci_method <- "bucketfill"
methods <- c("blp")
n_methods <- length(methods)

data_type <- "Scheetz2006"
alpha <- .2
lambda <- "cv"
ci_method <- "quantile"
args_list <- list(data = data_type, method = methods, lambda = lambda,
     ci_method = ci_method, nominal_coverage = alpha * 100)
params_grid <- expand.grid(args_list)

# Koussounadis2014, Scheetz2006, whoari
dat <- hdrm::readData("Scheetz2006")
# dat <- hdrm::readData("whoari")
dup <- duplicated(t(dat$X))
const <- apply(dat$X, 2, function(x) length(unique(x)) == 1)
dat$X <- dat$X[,!dup & !const]

if (any(!(methods %in% c("selectiveinference", "blp")))) {
  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min
  sigma2 <- cv_fit$cve[cv_fit$min]
}

cis <- list()
for (i in 1:n_methods) {
  print(methods[i])
  set.seed(my_seed)
  start <- Sys.time()
  if (methods[i] == "selectiveinference") {
    ### Selective Inference
    dat$X <- ncvreg::std(dat$X)
    res <- selective_inference(dat, estimate_sigma = FALSE)
  } else if (methods[i] == "blp") {
    ### HDI - Across a range of lambda values
    res <- blp(dat)
  } else {
    lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = TRUE, method = methods[i], nboot = nboot, lambda = lambda, sigma2 = sigma2)
    ci <- ci.boot.ncvreg(lassoboot, ci_method = ci_method) %>%
      mutate(method = methods[i])
    res <- list("confidence_interval" = ci, lambda, "sigma" = sqrt(sigma2))
  }
  print(as.numeric(difftime(Sys.time(), start, units = "secs")))

  save_objects(folder = rds_folder, args_list = args_list, overwrite = TRUE, res)
}

