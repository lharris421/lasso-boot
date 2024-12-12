## Setup
source("./scripts/setup.R")

libs <- c("ncvreg", "hdi", "dplyr", "hdrm")
evalFunc <- function() {
  devtools::load_all()
}
exports <- c("posterior", "pipe_ncvreg", "lp", "ci_full_cond", "soft_threshold", "firm_threshold_c")

methods <- methods[c("pipe_profile")]
simulation_info <- list(
  simulation_function = "gen_ortho",
  simulation_arguments = list(
    n = 100, p = 100,
    beta = c(-2, 2, -1, 1, -0.5, 0.5, -0.5, 0.5, rep(0, 92))
  )
)

## Specify number of iterations in run sim
run_sim(methods, simulation_info, parallel = TRUE, nClust = 5, libsEx = libs,
        clusterEv = evalFunc, clusterEx = exports, overwrite = TRUE)
