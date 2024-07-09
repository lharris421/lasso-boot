source("./scripts/setup/setup.R")
library(tictoc)
methods <- c("blp")
n_methods <- length(methods)

# data_type <- "brca1"
data_type <- "Scheetz2006"
# data_type <- "whoari"
alpha <- .2
lambda <- "cv"
args_list <- list(
  data = data_type,
  method = methods,
  lambda = lambda,
  nominal_coverage = (1-alpha) * 100,
  alpha = NA
)
params_grid <- expand.grid(args_list)

# Koussounadis2014, Scheetz2006, whoari
dat <- hdrm::readData("Scheetz2006")
dup <- duplicated(t(dat$X))
const <- apply(dat$X, 2, function(x) length(unique(x)) == 1)
dat$X <- dat$X[,!dup & !const]

set.seed(my_seed)
cv_fit <- cv.glmnet(dat$X, dat$y)
# cv.ncvreg(dat$X, dat$y, penalty ="lasso")
lambda <- cv_fit$lambda.min
sigma2 <- cv_fit$cvm[which.min(cv_fit$cvm)]

cis <- list()
for (i in 1:n_methods) {
  print(methods[i])
  set.seed(my_seed)
  start <- Sys.time()
  if (methods[i] == "selectiveinference") {
    ### Selective Inference
    dat$X <- ncvreg::std(dat$X)
    res <- selective_inference(dat, estimate_sigma = TRUE, lambda = lambda)
    print(res$lambda)
  } else if (methods[i] == "blp") {
    ### HDI - Across a range of lambda values
    # res <- blp(dat, boot.shortcut = TRUE, lambda = lambda)
    res <- blp(dat, boot.shortcut = TRUE)
    glmnet_fit <- glmnet(dat$X, dat$y)
    ests <- coef(glmnet_fit, s = lambda)[-1] ## put back to res$lambda if want selected by blp
    res$confidence_interval <- res$confidence_interval %>%
      mutate(estimate = ests)

  } else {
    lassoboot <- boot_ncvreg(dat$X, dat$y, verbose = TRUE, nboot = nboot, lambda = lambda, sigma2 = sigma2)
    ci <- ci.boot_ncvreg(lassoboot)
    res <- list("confidence_interval" = ci, lambda, "sigma" = sqrt(sigma2))
  }
  print(as.numeric(difftime(Sys.time(), start, units = "secs")))
  # args_list$method <- "blpmin"
  save_objects(folder = rds_path, args_list = args_list, overwrite = TRUE, res)
}

