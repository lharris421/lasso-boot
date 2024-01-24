source("./scripts/setup/setup.R")

# Koussounadis2014, Scheetz2006, whoari
# dat <- hdrm::readData(whoari)
dat <- hdrm::readData(Scheetz2006)
dup <- duplicated(t(dat$X))
const <- apply(dat$X, 2, function(x) length(unique(x)) == 1)
dat$X <- dat$X[,!dup & !const]
n <- nrow(dat$X)
p <- ncol(dat$X)

cv_fit <- cv.ncvreg(dat$X, dat$y, penalty = "lasso")
(lambda <- cv_fit$lambda.min)
(sigma2 <- cv_fit$cve[cv_fit$min])
lambdas <- rep(lambda, length(methods))
names(lambdas) <- methods

cis <- list()
for (i in 1:n_methods) {
  print(methods[i])
  lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = TRUE, quantiles = methods[i], nboot = nboot, lambda = lambda, sigma2 = sigma2)
  cis[[i]] <- ci.boot.ncvreg(lassoboot, method = method) %>%
    mutate(method = methods[i])
}

save(cis, lambdas, n, p, methods, file = glue("./rds/lassoboot_comparison_scheetz2006.rds"))
