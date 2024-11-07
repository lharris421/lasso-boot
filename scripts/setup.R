# Clear enviornment
rm(list=ls())

devtools::load_all()

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
library(tictoc)
library(dissertation)

methods <- list(
  "lasso_boot" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid")),
  "elastic_net" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid", enet_alpha = 0.8)),
  "traditional" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "traditional")),
  "posterior" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "posterior")),
  "selective_inference" = list(method = "selective_inference", method_arguments = list()),
  "lasso_proj_boot" = list(method = "blp", method_arguments = list()),
  "lasso_proj_boot_shortcut" = list(method = "blp", method_arguments = list(boot.shortcut = TRUE)),
  "ridge" = list(method = "ridge", method_arguments = list())
)
for (i in 1:length(methods)) {
  methods[[i]]$method_arguments["alpha"] <- 0.2
}
