## Setup
source("./scripts/setup.R")

# Ran w/o correlations between noise variables
# Ran w/ ar correlation = 0.5 between nopise variables

## Set up for all simulations
libs <- c("ncvreg", "hdi", "dplyr", "hdrm")
evalFunc <- function() {
  devtools::load_all()
}
exports <- c("posterior", "pipe_ncvreg", "lp", "ci_full_cond", "soft_threshold", "firm_threshold_c")

methods <- methods["pipe"]

simulation_info <- list(
  simulation_function = "gen_data_abn",
  simulation_arguments = list(
    n = 100, p = 100, a = 10, b = 4, rho = 0.5,
    beta = c(1 / (2^(0:9)) * c(1,-1),
             rep(0, 90))
    # rho.noise = 0.5, noise = "autoregressive"
  )
)

simulation_info <- list(
  simulation_function = "gen_data_abn",
  simulation_arguments = list(
    n = 100, p = 100, a = 10, b = 4, rho = 0.5,
    beta = c(rep(0.1, 10) * c(1,-1),
             rep(0, 90))
    # rho.noise = 0.5, noise = "autoregressive"
  )
)

run_sim(methods, simulation_info, parallel = TRUE, nClust = 5, libsEx = libs,
        clusterEv = evalFunc, clusterEx = exports)



