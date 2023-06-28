lp_dist <- function(b2, X, b1) {
  mean((X%*%b1 - X%*%b2[-1])^2)
}
min_lp_dist <- function(X, y, true_beta) {
  mod <- ncvreg(X, y, penalty = "lasso")
  err <- unlist(apply(mod$beta, 2, lp_dist, X, true_beta))
  return(mod$lambda[which.min(err)])
}
