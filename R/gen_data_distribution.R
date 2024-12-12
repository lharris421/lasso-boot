#' Title
#'
#' @param n
#' @param p
#' @param distribution
#' @param corr
#' @param rho
#' @param SNR
#'
#' @return
#' @export
#'
#' @examples
gen_data_distribution <- function(n = 100, p = 100, p1=floor(p/2),
                                  distribution = NULL, beta = NULL,
                                  corr = c("exchangeable", "autoregressive"),
                                  family = c("gaussian", "binomial", "poisson"),
                                  rho = 0, SNR = 1, sigma = 1, a = 1, b = 1,
                                  J = NULL, K = NULL, J1 = NULL, K1 = NULL, rho.g = NULL,
                                  rho.gz = NULL, ortho = FALSE) {

  corr <- match.arg(corr)
  family <- match.arg(family)

  if (!is.null(distribution) && distribution %in% c("whoari", "Scheetz2006")) {

    data <- do.call(hdrm::readData, list(name = distribution))
    dup <- duplicated(t(data$X))
    const <- apply(data$X, 2, function(x) length(unique(x)) == 1)
    data$X <- data$X[,!dup & !const]
    return(data)

  } else if (!is.null(distribution)) {
    if (distribution == "laplace") {
      betas <- rlaplace(p, rate = 1)
    } else if (distribution == "normal") {
      betas <- rnorm(p, sd = 2)
    } else if (distribution == "t") {
      betas <- rt(p, df = 4)
    } else if (distribution == "uniform") {
      betas <- runif(p, -1, 1)
    } else if (distribution == "beta") {
      betas <- rbeta(p, .1, .1) - .5
    } else if (distribution == "sparse 1") {
      betas <- c(rep(c(rep(0.5, 3), 1, 2), 2) * c(rep(1, 5), rep(-1, 5)), rep(0, 90))
    } else if (distribution == "sparse 2") {
      betas <- c(rnorm(30), rep(0, 70))
    } else if (distribution == "sparse 3") {
      betas <- c(rnorm(50), rep(0, 50))
    } else if (distribution == "sparse 4") {
      betas <- c(3, rep(0, 99))
    } else if (distribution == "epsilon_conundrum") {
      betas <- c(-2, -1, -0.5, 0.5, 1, 2, 0 + 1e-9 * sample(c(-1, 1), 94, replace = TRUE))
    } else if (distribution == "high_corr") {
      betas <- 1
    } else if (distribution %in% c("abn", "group")) {
      betas <- beta
    } else if (distribution == "ldpe") {
      which_nonzero <- c(1500, 1800, 2100, 2400, 2700, 3000)
      true_lambda <- sqrt((2/n)*log(p))
      betas <- 3*true_lambda / ((1:p)^a)
      betas[which_nonzero] <- 3*true_lambda
    }
  }


  ## Generate data
  if (is.null(beta)) {
    if (!is.null(distribution) && distribution == "high_corr") {
      data <- gen_data_abn(
        n = n, p = 100, a = 1, b = 1,
        beta = betas, rho = 0.99
      )
    } else if (!is.null(distribution) && distribution == "ldpe") {
      data <- gen_data(n, p, beta = betas, corr = corr, rho = rho)
    } else if (is.null(distribution)) {
      data <- gen_data_snr(
        n = n, p = p, p1 = p1, family = family,
        corr = corr, rho = rho, SNR = SNR, sigma = sigma
      )
    } else {
      data <- gen_data_snr(
        n = n, p = p, p1 = p1, family = family,
        beta = betas, corr = corr, rho = rho, SNR = SNR, sigma = sigma
      )
    }
  } else {
    if (!is.null(distribution) && distribution == "abn") {
      data <- gen_data_abn(
        n = n, p = p, a = a, b = b,
        beta = betas, rho = rho
      )
    } else if (!is.null(distribution) && distribution == "group") {
        data <- gen_data_grp(
          n = n, J = J, K = K, beta = betas,
          J1 = J1, K1 = K1, rho.g = rho.g, rho.gz = rho.gz
        )
    } else if (ortho) {
      data <- gen_ortho(
        n = n, p = length(beta),
        beta = beta
      )
    } else {
      data <- gen_data(
        n = n, p = length(beta), family = family,
        beta = beta, corr = corr, rho = rho
      )
    }
  }

  return(data)

}

