unloadNamespace("hdrm")
unloadNamespace("ncvreg")
rm(list=ls())
.libPaths("./local")
library(ncvreg)
library(hdrm)
library(glue)

library(glmnet)
library(selectiveInference)
library(hdi)

library(dplyr)
library(stringr)
library(tidyr)
library(monomvn)

devtools::load_all()

my_seed <- 189807771
set.seed(my_seed)

nboot <- 1000

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}

ci_method <- "quantile"
methods <- c("traditional", "sample", "debiased", "acceptreject", "zerosample1", "zerosample2")
# methods <- c("traditional", "sample", "debiased", "zerosample2")
n_methods <- length(methods)

## Alt gen data
gen_data_snr <- function(n, p, p1=floor(p/2), beta, family=c("gaussian","binomial"), SNR=1,
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
  # sigma <- ifelse(missing(beta) || length(beta)==1, 1, sqrt(var(X %*% beta)/SNR))
  sigma <- ifelse(missing(beta) || length(beta)==1, 1, sqrt(drop(crossprod(beta))/SNR))
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

gen_y <- function(eta, family=c("gaussian", "binomial"), sigma=1) {
  family=match.arg(family)
  n <- length(eta)
  if (family=="gaussian") {
    rnorm(n, mean=eta, sd=sigma)
  } else if (family=="binomial") {
    pi. <- exp(eta)/(1+exp(eta))
    pi.[eta > log(.9999/.0001)] <- 1
    pi.[eta < log(.0001/.9999)] <- 0
    rbinom(n,1,pi.)
  }
}

gen_data_abn2 <- function(n=100, p=60, a=6, b=2, rho=0.5, family=c("gaussian", "binomial"), signal=c('homogeneous', 'heterogeneous'), noise=c('exchangeable', 'autoregressive'),
                         rho.noise=0, beta, SNR=1) {

  family <- match.arg(family)
  signal <- match.arg(signal)
  noise <- match.arg(noise)
  K <- b + 1

  # Gen X
  columns <- list()
  for (i in 1:a) {
    columns[[i]] <- gen_x(n, K, rho, noise)
  }
  columns[[a + 1]] <- gen_x(n, p - a*K, rho.noise, noise)
  X <- do.call(cbind, columns)

  # Gen beta
  if (missing(beta) || length(beta)==1) {
    bb <- c(-1,1)[(1:a)%%2+1]
    if (missing(beta)) {
      if (signal=="heterogeneous") bb <- bb*(a:1)
      bbb <- numeric(p)
      bbb[((1:a)-1)*K+1] <- bb
      beta <- bbb*sqrt(SNR)/sqrt(drop(crossprod(bbb)))
    } else {
      bbb <- numeric(p)
      bbb[((1:a)-1)*K+1] <- bb
      beta <- bbb*beta
    }
  } else {
    bb <- beta
    beta <- numeric(p)
    beta[((1:a)-1)*K+1] <- bb
  }

  beta <- bbb
  print(bbb)
  print(drop(crossprod(beta)))
  # Gen y
  y <- gen_y(X%*%beta, family=family, sigma=sqrt(drop(crossprod(beta))/SNR))


  # Return
  varType <- vector("character", p)
  varType[((1:a)-1)*K+1] <- "A"
  for (j in 1:b) {
    varType[((1:a)-1)*K+1+j] <- "B"
  }
  varType[(a*K+1):p] <- "N"
  varLab <- vector("character", p)
  varLab[varType=="A"] <- paste0("A", 1:sum(varType=="A"))
  varLab[varType=="B"] <- paste0("B", 1:sum(varType=="B"))
  varLab[varType=="N"] <- paste0("N", 1:sum(varType=="N"))
  colnames(X) <- varLab
  names(beta) <- varLab
  list(X=X, y=y, beta=beta, family=family, varType=varType)
}


