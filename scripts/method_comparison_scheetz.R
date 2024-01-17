rm(list=ls())
unloadNamespace("ncvreg")
.libPaths("./local")
library(ncvreg)
library(ggplot2)
library(dplyr)

my_seed <- 189807771
set.seed(my_seed)

plot_res <- list()
# Koussounadis2014, Scheetz2006
dat <- hdrm::readData(Scheetz2006)
n <- nrow(dat$X)
p <- ncol(dat$X)

### Lasso-boot
lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, nboot = 1000, quantiles = "disturbed")
lassoboot_ci <- ci.boot.ncvreg(lassoboot)
lassoboot_lam <- lassoboot$lamdba

plot_res <- list(lassoboot_ci, lassoboot_lam, n, p)

plot(lassoboot, n = 50)
plot(lassoboot, n = 500)

save(plot_res, file = "./rds/method_comparison_Scheetz2006.rds")

diffs <- numeric(p)
for (i in 1:p) {
  ex <- lassoboot$draws[,i]

  rsa <- create_rsa(ex[1:100], quantile = .9)
  rsa <- update(rsa, ex[101:1000])
  diffs[i] <- quantile(ex, .9) - rsa$xi
}

mean(diffs^2)
mean(abs(diffs))
