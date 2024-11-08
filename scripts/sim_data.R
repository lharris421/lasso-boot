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
#   simulation_function = "gen_data",
#   simulation_arguments = list(
#     n = 50, p = 100,
#     beta = c(rep(0.5, 50), rep(0, 50))
#     # corr = "autoregressive", rho = 0.5
#   )
# )

simulation_info <- list(
  simulation_function = "gen_data",
  simulation_arguments = list(
    n = 100, p = 100,
    beta = c(-2, 2, -1, 1, -0.5, 0.5, -0.5, 0.5, rep(0, 92))
  )
)

# set.seed(07803715)
# exp_beta <- rexp(50, 1) * sample(c(-1, 1), 50, replace = TRUE)
# simulation_info <- list(
#   simulation_function = "gen_data",
#   simulation_arguments = list(
#     n = 40, p = 100,
#     beta = c(exp_beta, rep(0, 50))
#   )
# )

run_sim(methods, simulation_info, parallel = TRUE, nClust = 5, libsEx = libs,
        clusterEv = evalFunc, clusterEx = exports)

