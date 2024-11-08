## Setup
source("./scripts/setup.R")

## Set up for all simulations
libs <- c("ncvreg", "hdi", "dplyr", "hdrm")
evalFunc <- function() {
  devtools::load_all()
}
exports <- c("posterior", "pipe_ncvreg", "lp", "ci_full_cond", "soft_threshold",
             "firm_threshold_c", "gen_data_ldpe")

methods <- methods[c("mcp_cond")]
simulation_info <- list(
  simulation_function = "gen_data_ldpe",
  simulation_arguments = list(
    n = 200, p = 3000,
    corr = "autoregressive",
    rho = 0.8, a = 1
  )
)

run_sim(methods, simulation_info, parallel = TRUE, nClust = 3, libsEx = libs,
        clusterEv = evalFunc, clusterEx = exports)

