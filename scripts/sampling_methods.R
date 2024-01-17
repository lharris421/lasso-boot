## Sampling methods
set.seed(189807771)
n <- 1000

## paramter of interest
mn <- rnorm(1e2, 0, .1)

## Original sample
os <- rnorm(1e2, mn, 1)

## Simple
draws1 <- numeric(n)
for (i in 1:n) {
  samp <- sample(os, size = length(os), replace = TRUE)
  ps <- runif(1)
  draws1[i] <- qnorm(ps, mean = mean(samp), sd = sd(samp) / sqrt(length(samp)))
}
dens1 <- density(draws1, bw = .1)

## Traditional
draws2 <- numeric(n)
for (i in 1:n) {
  draws2[i] <- mean(sample(os, size = length(os), replace = TRUE))
}
dens2 <- density(draws2, bw = .1)

## Accept / reject -> simplify
draws3 <- numeric(n)
for (i in 1:n) {
  samp <- sample(os, size = length(os), replace = TRUE)
  draw <- runif(1, -3, 3) * (sd(samp) / sqrt(length(samp))) + mean(samp)
  thresh <- runif(1)
  all_dens <- dnorm(draw, mean = mean(samp), sd = sd(samp)/ sqrt(length(samp))) / dnorm(mean(samp), mean(samp), sd(samp)/ sqrt(length(samp)))
  while(any(all_dens < thresh)) {
    failed <- draw[all_dens < thresh]
    draw[all_dens < thresh] <- runif(length(failed), -abs(mean(samp) - failed), abs(mean(samp) - failed)) * (sd(samp) / sqrt(length(samp))) + mean(samp)
    thresh[all_dens < thresh] <- runif(length(failed))
    all_dens <- dnorm(draw, mean = mean(samp), sd = sd(samp)/ sqrt(length(samp))) / dnorm(mean(samp), mean(samp), sd(samp)/ sqrt(length(samp)))
  }
  draws3[i] <- draw
}
dens3 <- density(draws3, bw = .1)

plot(x = dens1$x, y = dens1$y / max(dens1$y), type = "l", col = "red")
lines(x = dens2$x, y = dens2$y / max(dens2$y), col = "blue")
lines(x = dens3$x, y = dens3$y / max(dens3$y), col = "green")
lines(x = seq(-4, 4, by = 0.01), y = dnorm(seq(-4, 4, by = 0.01), sd = .1) / dnorm(0, sd = .1))
