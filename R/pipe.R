pipe_ncvreg <- function(
    X, y, fit, lambda, sigma,
    family = "gaussian",
    penalty = "lasso",
    gamma = switch(penalty, SCAD = 3.7, 3),
    enet_alpha = 1, alpha = 0.05,
    posterior = FALSE
) {

  res <- ncvreg::pipe(
    X = X, y = y, fit = fit, lambda = lambda, sigma = sigma,
    family = family, penalty = penalty, gamma = gamma,
    alpha = enet_alpha, level = 1 - alpha,
    posterior = posterior
  )

  return(res)

}
