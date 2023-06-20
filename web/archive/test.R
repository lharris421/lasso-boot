## Mixture attempt 1
q1 <- qunif((1:3000)/3001)
q2 <- qnorm((1:7000)/7001, mean = 3, sd = 2)
quantile(c(q1, q2), c(.025, .975))

r1 <- runif(1000000)
r2 <- rnorm(1000000, 3, 2)
quantile(sample(c(r1, r2), size = 1000000, prob = c(rep(.3, 1000000), rep(.7, 1000000)), replace = TRUE), c(.025, .975))


