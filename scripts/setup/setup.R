# Clear enviornment
rm(list=ls())

# Load packages
library(hdrm)
library(glue)
library(glmnet)
library(selectiveInference)
library(hdi)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(indexr)
library(digest)

devtools::load_all()

ci_method <- "quantile"
methods <- c("traditional", "sample", "debiased", "zerosample2")
n_methods <- length(methods)

nboot <- 1000

my_seed <- 189807771
set.seed(my_seed)
