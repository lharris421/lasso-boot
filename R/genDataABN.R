#' Title
#'
#' @param n holder
#' @param p holder
#' @param a holder
#' @param b holder
#' @param rho holder
#' @param family holder
#' @param signal holder
#' @param noise holder
#' @param rho.noise holder
#' @param beta holder
#' @param SNR holder
#'
#' @return
#' @export
#'
#' @examples
genDataABN <-
  function (n = 100,
            p = 60,
            a = 6,
            b = 2,
            rho = 0.5,
            family = c("gaussian",
                       "binomial"),
            signal = c("homogeneous", "heterogeneous"),
            noise = c("exchangeable", "autoregressive"),
            rho.noise = 0,
            beta,
            SNR = 1, sgm = 1)
  {
    family <- match.arg(family)
    noise <- match.arg(noise)
    signal <- match.arg(signal)
    K <- b + 1
    sigmaList <- vector("list", a + 1)
    for (i in 1:a) {
      sigmaList[[i]] <- matrix(rho, K, K) + (1 - rho) * diag(K)
    }
    sigmaList[[a + 1]] <- hdrm:::genS(p - K * a, rho.noise, noise)
    S <- Matrix::.bdiag(sigmaList)
    X <- hdrm:::genX(n, p, S)
    if (missing(beta) || length(beta) == 1) {
      bb <- c(-1, 1)[(1:a) %% 2 + 1]
      if (missing(beta)) {
        if (signal == "heterogeneous")
          bb <- bb * (a:1)
        bbb <- numeric(p)
        bbb[((1:a) - 1) * K + 1] <- bb
        beta <-
          bbb * sqrt(SNR) / sqrt(Matrix::drop(Matrix::crossprod(bbb, S) %*% bbb))
      }
      else {
        bbb <- numeric(p)
        bbb[((1:a) - 1) * K + 1] <- bb
        beta <- bbb * beta
      }
    }
    else {
      bb <- beta
      beta <- numeric(p)
      beta[((1:a) - 1) * K + 1] <- bb
    }
    y <- hdrm:::genY(X %*% beta, family = family, sigma = sgm)
    varType <- vector("character", p)
    varType[((1:a) - 1) * K + 1] <- "A"
    if (b > 0) {
      for (j in 1:b) {
        varType[((1:a) - 1) * K + 1 + j] <- "B"
      }
    }
    if (a + b < p)
      varType[(a * K + 1):p] <- "N"
    varLab <- vector("character", p)
    varLab[varType == "A"] <- paste0("A", 1:sum(varType == "A"))
    varLab[varType == "B"] <- paste0("B", 1:sum(varType == "B"))
    varLab[varType == "N"] <- paste0("N", 1:sum(varType == "N"))
    colnames(X) <- varLab
    names(beta) <- varLab
    list(
      X = X,
      y = y,
      beta = beta,
      family = family,
      varType = varType
    )
  }
