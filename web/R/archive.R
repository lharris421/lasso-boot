## Alternative ways of calculating lambda

### Estimation Consistency
prob <- .99
p <- 100
n <- 160
c <- log((prob - 1) / -2) * (-2 / log(p)) + 2
(lambda <- 2*sqrt(c*log(p)/n))

eval_scenario("ec", type1 = "pairs", type2 = "percentile", n = n, p = p, a = a, b = b, betas = betas, lambda = lambda)
