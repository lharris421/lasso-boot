gen_data_sigma <- function(n, p, p1=floor(p/2), beta, family=c("gaussian","binomial","poisson"), SNR=1, sigma = 1,
                     signal = c("homogeneous","heterogeneous"), corr=c("exchangeable", "autoregressive"),
                     rho = 0) {
  family <- match.arg(family)
  signal <- match.arg(signal)
  corr <- match.arg(corr)

  # Gen X
  X <- gen_x(n, p, rho, corr)

  # Gen beta
  if (missing(beta) || length(beta)==1) {
    j <- 1:p
    s <- c(-1,1)[j%%2+1]
    b <- (j <= p1)
    if (missing(beta)) {
      if (signal=="heterogeneous") b <- b*rev(j)
      b <- b*s
      beta <- b*sqrt(SNR)/sqrt(drop(crossprod(b)))
      #beta <- b*sqrt(SNR)/sqrt(calc_bsb(b[1:p1], rho, corr))
    } else {
      beta <- b*s*beta
    }
  }

  # Gen y
  y <- gen_y(X%*%beta, family=family, sigma=sigma)

  # Label and return
  w <- 1 + floor(log10(p))
  vlab <- paste0('V', formatC(1:p, format='d', width=w, flag='0'))
  colnames(X) <- names(beta) <- vlab
  list(X=X, y=y, beta=beta, family=family)
}

gen_x <- function(n, p, rho, corr=c('exchangeable', 'autoregressive')) {
  corr <- match.arg(corr)
  if (corr == 'exchangeable') {
    z <- rnorm(n)
    sqrt(rho)*z + sqrt(1-rho) * matrix(rnorm(n*p), n, p)
  } else if (corr == 'autoregressive') {
    Z <- cbind(rnorm(n), matrix(rnorm(n*(p-1), sd=sqrt(1-rho^2)), n, p-1))
    apply(Z, 1, stats::filter, filter=rho, method='recursive') |> t()
  }
}

calc_bsb <- function(b, rho, corr) {
  if (corr == 'exchangeable') {
    sum(rho*tcrossprod(b)) + (1-rho)*crossprod(b) |> drop()
  } else if (corr == 'autoregressive') {
    out <- crossprod(b)
    bb <- tcrossprod(b)
    for (j in 1:min(10, length(b)-1)) {
      out <- out + 2 * rho^j * sum(Matrix::band(bb, j, j))
    }
    drop(out)
  }
}

gen_y <- function(eta, family = c("gaussian", "binomial", "poisson"), sigma = 1) {
  family <- match.arg(family)
  n <- length(eta)
  if (family == "gaussian") {
    rnorm(n, mean = eta, sd = sigma)
  } else if (family == "binomial") {
    pi. <- 1 / (1 + exp(-eta))
    pi.[eta > log(.9999 / .0001)] <- 1
    pi.[eta < log(.0001 / .9999)] <- 0
    rbinom(n, 1, pi.)
  } else if (family == "poisson") {
    mu <- exp(eta)
    rpois(n, lambda = mu)
  }
}
