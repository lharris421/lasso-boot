set.seed(my_seed)
niter <- 100
coverage <- numeric(niter)
n <- 1000
p <- 100
alpha <- .2
data_type <- "normal"
method <- "lm"
ci_method <- "quantile"
SNR <- 1
corr <- "autoregressive"
rho <- 0.8
sd <- 1
res <- list()
for (j in 1:niter) {
  print(j)
  dat <- gen_data_snr(n = n, p = p, p1 = p, beta = rnorm(p, mean = 0, sd = sd), corr = corr, rho = rho)
  df <- data.frame(dat$y, dat$X)
  colnames(df) <- c("y", names(dat$beta))

  draws <- matrix(nrow = nboot, ncol = p)
  for (i in 1:nboot) {
    idx_new <- sample(1:nrow(df), replace = TRUE)
    fit <- lm(y ~ ., data = df[idx_new,])
    draws[i,] <- coef(fit)[-1]
  }
  failed <- (max(apply(draws, 2, function(x) sum(is.na(x)))) / nboot)
  print(glue::glue("Highest failure percentage: {failed*100}%"))
  if(failed > .1) stop("Too unstable, increase n")
  lowers <- apply(draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(draws, 2, function(x) quantile(x, 1 - (alpha / 2), na.rm = TRUE))

  coverage[j] <- mean(lowers <= dat$beta & uppers >= dat$beta)
  print(glue::glue("Coverage for simulation {j}: {coverage[j]}"))

  res[[j]] <- data.frame(estimate = coef(lm(y ~., df))[-1], variable = names(dat$beta), lower = lowers, upper = uppers, ci_method = "quantile", truth = dat$beta)
}

folder <- "./rds/"
summary(coverage)
args_list <- list(data = data_type,
                  n = n,
                  snr = SNR,
                  sd = ifelse(data_type == "normal", sd, NA),
                  rate = ifelse(data_type == "laplace", rt, NA),
                  a = ifelse(data_type == "abn", a, NA),
                  b = ifelse(data_type == "abn", a, NA),
                  correlation_structure = corr,
                  correlation = rho,
                  correlation_noise = ifelse(data_type == "abn", rho.noise, NA),
                  method = method,
                  ci_method = ci_method,
                  nominal_coverage = alpha * 100,
                  modifier = NA,
                  lambda = "cv",
                  p = p)
save_objects(folder = folder, args_list = args_list, res, save_method = "rda")
