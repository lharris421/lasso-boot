source("./scripts/setup/setup.R")

quantiles <- "disturbed"
method <- "quantile"

plot_res <- list()

dat <- hdrm::readData(whoari)
dup <- duplicated(t(dat$X))
const <- apply(dat$X, 2, function(x) length(unique(x)) == 1)
dat$X <- dat$X[,!dup & !const]
dat$X <- ncvreg::std(dat$X)
# dat$y <- ncvreg::std(dat$y)
n <- nrow(dat$X)
p <- ncol(dat$X)

### Selective Inference
cv_res <- cv.glmnet(dat$X, dat$y, standardize = FALSE)
si_lam <- lam <- cv_res$lambda.min

fit <- cv_res$glmnet.fit
b <- coef(fit, s = lam, exact = TRUE)[-1]
names(b) <- colnames(dat$X)
tryCatch({
  # sh <- estimateSigma(dat$X, dat$y)$sigmahat
  # print(sh)
  res <- fixedLassoInf(dat$X, dat$y, b, lam*length(dat$y), alpha = .2)
  B <- res$ci
  rownames(B) <- names(res$vars)
  colnames(B) <- c("lower", "upper")
  print(B)
  # B <- B[is.finite(B[,2]) & is.finite(B[,3]),-4]
  si_ci <- B %>%
    data.frame(method = "Selective Inference", variable = rownames(B)) %>%
    mutate(estimate = b[rownames(B)])
}, error = function(e) {
  si_ci <- NA
  print(e)
}
)

### HDI - Across a range of lambda values
fit.lasso.allinfo <- boot.lasso.proj(dat$X, dat$y, return.bootdist = TRUE, B = nboot, boot.shortcut = TRUE)
ci_hdi <- confint(fit.lasso.allinfo, level = 0.8)

hdi_ci <- ci_hdi %>%
  data.frame(method = "BLP", variable = rownames(ci_hdi)) %>%
  mutate(estimate = fit.lasso.allinfo$bhat)
hdi_lam <- fit.lasso.allinfo$lambda

### Lasso-boot
lassoboot <- boot.ncvreg(dat$X, dat$y, verbose = FALSE, quantiles = quantiles, nboot = nboot)
lassoboot_ci <- ci.boot.ncvreg(lassoboot, method = method, original_data = dat)
lassoboot_lam <- lassoboot$lamdba

plot_res <- list(si_ci, hdi_ci, lassoboot_ci, si_lam, hdi_lam, lassoboot_lam, n, p)

save(plot_res, file = glue("./rds/method_comparison_whoari_{quantiles}_{method}_span1.rds"))
