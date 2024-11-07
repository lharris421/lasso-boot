## Setup
source("./scripts/setup.R")

## Set up for all simulations
libs <- c("ncvreg", "hdi", "dplyr", "hdrm", "dissertation", "selectiveInference", "glmnet")
evalFunc <- function() {
  devtools::load_all()
}
# exports <- c("boot_ncv", "lp", "blp", "selective_inference")


methods <- methods[c("lasso_boot", "lasso_proj_boot")]
# methods <- methods[c("selective_inference")]
simulation_info <- list(
  simulation_function = "gen_data_distribution",
  simulation_arguments = list(
    n = 100, p = 100, SNR = 1,
    distribution = "laplace"
  ),
  script_name ="distributions"
)

run_sim(methods, simulation_info, parallel = TRUE, nClust = 2, libsEx = libs,
        clusterEv = evalFunc, mimic_script = TRUE)
