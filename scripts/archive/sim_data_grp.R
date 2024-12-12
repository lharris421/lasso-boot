## Setup
source("./scripts/setup.R")

## Set up for all simulations
libs <- c("ncvreg", "hdi", "dplyr", "hdrm")
evalFunc <- function() {
  devtools::load_all()
}
exports <- c("posterior", "pipe_ncvreg", "lp", "ci_full_cond", "soft_threshold", "firm_threshold_c")

methods <- methods[c("pipe_profile")]
# simulation_info <- list(
#   simulation_function = "gen_data_grp",
#   simulation_arguments = list(
#     n = 100, J = 25, K = 4, beta = c(1, 1, 1, rep(0, 97)),
#     J1 = 1, K1 = 3, rho.g = 0.5, rho.gz = 0 # K1 doesn't actually matter here (not used)
#   )
# )

simulation_info <- list(
  simulation_function = "gen_data_grp",
  simulation_arguments = list(
    n = 100, J = 20, K = 5, beta = 0.5,
    J1 = 2, K1 = 3, SNR = 1, rho.g = 0.5
  )
)


# simulation_info <- list(
#   simulation_function = "gen_data_grp",
#   simulation_arguments = list(
#     n = 100, J = 25, K = 4, beta = 0.5,
#     J1 = 1, K1 = 4, rho.g = 0.5
#   )
# )

run_sim(methods, simulation_info, parallel = TRUE, nClust = 5, libsEx = libs,
        clusterEv = evalFunc, clusterEx = exports)

