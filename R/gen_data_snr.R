#' Title
#'
#' @param n
#' @param p
#' @param p1
#' @param beta
#' @param family
#' @param SNR
#' @param signal
#' @param corr
#' @param rho
#'
#' @return
#' @export
#'
#' @examples
gen_data_snr <- function(n, p, p1=floor(p/2), beta, family=c("gaussian","binomial","poisson"), SNR=1, sigma = 1,
                         signal = c("homogeneous","heterogeneous"), corr=c("exchangeable", "autoregressive"),
                         rho = 0) {

  family <- match.arg(family)
  signal <- match.arg(signal)
  corr <- match.arg(corr)

  if (!missing(beta)) beta <- (beta / sqrt(drop(crossprod(beta)))) * sqrt(SNR) * sigma

  if (!missing(beta)) {
    gen_data_sigma(n = n, p = p, p1 = p1, beta = beta, family = family, SNR = SNR, sigma = sigma,
                   signal = signal, corr = corr, rho = rho)
  } else {
    gen_data_sigma(n = n, p = p, p1 = p1, family = family, SNR = SNR, sigma = sigma,
                   signal = signal, corr = corr, rho = rho)
  }

}
