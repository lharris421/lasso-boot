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

method <- "quantile"
methods <- c("traditional", "sample", "debiased", "acceptreject", "zerosample1", "zerosample2")
# methods <- c("traditional", "sample", "debiased", "zerosample2")
n_methods <- length(methods)
