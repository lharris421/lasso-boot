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
gen_data_distribution <- function(n = 100, p = 100,
                                  distribution = NULL, beta = NULL,
                                  corr = c("exchangeable", "autoregressive"),
                                  rho = 0, SNR = 1) {

  corr <- match.arg(corr)

  if (!is.null(distribution)) {

    if (distribution %in% c("whoari", "Scheetz2006")) {
      data <- do.call(hdrm::readData, list(name = distribution))
      dup <- duplicated(t(data$X))
      const <- apply(data$X, 2, function(x) length(unique(x)) == 1)
      data$X <- data$X[,!dup & !const]
      return(data)
    } else {
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
      }
    }
  } else {
    betas <- beta
  }

  if (is.null(beta)) {
    if (distribution %in% c("high_corr", "abn")) {
      data <- gen_data_abn(
        n = n, p = 100, a = 1, b = 1,
        beta = betas, rho = 0.99, SNR = SNR
      )
      return(data)
    } else {
      data <- gen_data_snr(
        n = n, p = length(betas),
        beta = betas, corr = corr, rho = rho, SNR = SNR
      )
      return(data)
    }
  } else {
    data <- gen_data(
      n = n, p = length(betas),
      beta = betas, corr = corr, rho = rho
    )
    return(data)
  }

}

