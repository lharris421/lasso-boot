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
library(ncvreg)
library(ggplot2)
library(magrittr)
library(parallel)

methods <- list(
  "lasso_boot" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid")),
  "mcp_boot" = list(method = "boot_ncv", method_arguments = list(penalty = "MCP", submethod = "hybrid")),
  "elastic_net" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid", enet_alpha = 0.8)),
  "traditional" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "traditional")),
  "posterior" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "posterior")),
  "debiased" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "debiased")),
  "debiased_fma" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "debiased", reselect_lambda = TRUE)),
  "selective_inference" = list(method = "selective_inference", method_arguments = list()),
  "lasso_proj_boot" = list(method = "blp", method_arguments = list()),
  "lasso_proj_boot_shortcut" = list(method = "blp", method_arguments = list(boot.shortcut = TRUE)),
  "ridge" = list(method = "ridge", method_arguments = list()),
  "lasso" = list(method = "posterior", method_arguments = list(penalty = "lasso")),
  "mcp"   = list(method = "posterior", method_arguments = list(penalty = "MCP")),
  "lasso_relaxed" = list(method = "posterior", method_arguments = list(penalty = "lasso", relaxed = TRUE)),
  "mcp_relaxed"   = list(method = "posterior", method_arguments = list(penalty = "MCP", relaxed = TRUE)),
  "normal_approx"  = list(method = "pipe_ncvreg", method_arguments = list(original_n = TRUE)),
  "pipe"  = list(method = "pipe_ncvreg", method_arguments = list()),
  "pipe_relaxed"  = list(method = "pipe_ncvreg", method_arguments = list(relaxed = TRUE)),
  "pipe_mcp"  = list(method = "pipe_ncvreg", method_arguments = list(penalty = "MCP"))
)
for (i in 1:length(methods)) {
  methods[[i]]$method_arguments["alpha"] <- 0.2
}
