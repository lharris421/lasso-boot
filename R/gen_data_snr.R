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
gen_data_snr <- function(n, p, p1=floor(p/2), beta, family=c("gaussian","binomial","hetero"), SNR=1,
                         signal = c("homogeneous","heterogeneous"), corr=c("exchangeable", "autoregressive"),
                         rho = 0) {

  family <- match.arg(family)
  signal <- match.arg(signal)
  corr <- match.arg(corr)

  beta <- (beta / sqrt(drop(crossprod(beta)))) * sqrt(SNR)

  gen_data(n = n, p = p, p1 = p1, beta = beta, family = family, SNR = SNR,
           signal = signal, corr = corr, rho = rho)

}
