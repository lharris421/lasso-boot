source("./scripts/setup/setup.R")

# ci_method <- "bucketfill"
methods <- c("zerosample1", "zerosample2", "sample", "debiased", "acceptreject", "traditional", "blp", "selective_inference")
n_methods <- length(methods)

# dataset <- "Scheetz2006"
dataset <- "whoari"
# Koussounadis2014, Scheetz2006, whoari
# dat <- hdrm::readData("Scheetz2006")
dat <- hdrm::readData("whoari")
dup <- duplicated(t(dat$X))
const <- apply(dat$X, 2, function(x) length(unique(x)) == 1)
dat$X <- dat$X[,!dup & !const]
dat$X <- ncvreg::std(dat$X)

if (!(any(methods) %in% c("selective_inference", "blp"))) {
  cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
  lambda <- cv_fit$lambda.min
  sigma2 <- cv_fit$cve[cv_fit$min]
}

cis <- list()
for (i in 1:n_methods) {
  print(methods[i])
  set.seed(my_seed)
  start <- Sys.time()
  if (methods[i] == "selective_inference") {
    ### Selective Inference
    res <- selective_inference(dat, estimate_sigma = TRUE)
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

  save(res, file = glue("./rds/{dataset}_{methods[i]}.rds"))
}

