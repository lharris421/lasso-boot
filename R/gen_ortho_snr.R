#' Title
#'
#' @param n
#' @param p
#' @param p1
#' @param beta
#' @param family
#' @param SNR
#' @param signal
#'
#' @return
#' @export
#'
#' @examples
gen_ortho_snr <- function(n, p, p1=floor(p/2), beta, family=c("gaussian","binomial"), SNR=1,
                          signal = c("homogeneous","heterogeneous")) {

  beta <- beta*sqrt(SNR)/sqrt(drop(crossprod(beta)))

  gen_ortho(n = n, p = p, p1 = p1, beta = beta, family = family, SNR = SNR,
            signal - signal)

}
