#unloadNamespace("ncvreg")
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

my_seed <- 189807771
set.seed(my_seed)

nboot <- 1000

rlaplace <- function(n, rate = 1) {
  rexp(n, rate) * sample(c(-1, 1), n, replace = TRUE)
}
