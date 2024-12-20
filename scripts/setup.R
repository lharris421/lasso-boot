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
  "lasso_boot_reed" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid", sigma2_reed = TRUE)),
  "mcp_boot" = list(method = "boot_ncv", method_arguments = list(penalty = "MCP", submethod = "hybrid")),
  "elastic_net" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "hybrid", enet_alpha = 0.8)),
  "traditional" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "traditional")),
  "posterior" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "posterior")),
  "debiased" = list(method = "boot_ncv", method_arguments = list(penalty = "lasso", submethod = "debiased")),
  "selective_inference" = list(method = "selective_inference", method_arguments = list()),
  "lasso_proj_boot" = list(method = "blp", method_arguments = list()),
  "lasso_proj_boot_shortcut" = list(method = "blp", method_arguments = list(boot.shortcut = TRUE)),
  "lasso_proj" = list(method = "lp", method_arguments = list()),
  "ridge" = list(method = "ridge", method_arguments = list()),
  # "lasso" = list(method = "posterior", method_arguments = list(penalty = "lasso")),
  "lasso"  = list(method = "pipe_ncvreg", method_arguments = list(posterior = TRUE)),
  "mcp"   = list(method = "posterior", method_arguments = list(penalty = "MCP")),
  "lasso_relaxed" = list(method = "posterior", method_arguments = list(penalty = "lasso", relaxed = TRUE)),
  "mcp_relaxed"   = list(method = "posterior", method_arguments = list(penalty = "MCP", relaxed = TRUE)),
  "normal_approx"  = list(method = "pipe_ncvreg", method_arguments = list(original_n = TRUE)),
  # "pipe"  = list(method = "pipe", method_arguments = list(level = 0.95, penalty = "lasso")),
  "pipe"  = list(method = "pipe_ncvreg", method_arguments = list()),
  "relaxed_lasso"  = list(method = "pipe_ncvreg", method_arguments = list(relaxed = TRUE)),
  "pipe_mcp"  = list(method = "pipe_ncvreg", method_arguments = list(penalty = "MCP"))
  # "pipe_poisson" = list(method = "pipe", method_arguments = list(family = "poisson", level = 0.95)),
  # "pipe_binomial" = list(method = "pipe", method_arguments = list(family = "binomial", level = 0.95))
)
