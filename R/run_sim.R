#' Title
#'
#' @param methods
#' @param simulation_info
#' @param iterations
#' @param overwrite
#' @param seed
#' @param parallel
#' @param nClust
#' @param libsEx
#' @param clusterEv
#' @param clusterEx
#'
#' @return
#' @export
#'
#' @examples
run_sim <- function(methods, simulation_info, iterations = 1000, overwrite = TRUE, seed = 1234, parallel = FALSE,
                    nClust = 2, libsEx = NULL, clusterEv = NULL, clusterEx = NULL, mimic_script = FALSE) {

  # Combine the method information and simulation information
  simulation_arguments <- list(
    "methods" = methods,
    "simulation_function" = simulation_info$simulation_function,
    "simulation_arguments" = simulation_info$simulation_arguments,
    "iterations" = iterations,
    "seed" = seed
  )

  if (mimic_script) {
    sim_args <- c(simulation_info, "seed" = seed, "iterations" = iterations)
  }

  # Run the simulation
  if (parallel) {
    simulation_arguments <- c(
      simulation_arguments,
      list("nClust" = 2, "libsEx" = libsEx, "clusterEv" = clusterEv, "clusterEx" = clusterEx)
    )

    simulation_results <- do.call(sim_parallel, simulation_arguments)

  } else {
    simulation_results <- do.call(sim, simulation_arguments)
  }

  if (!mimic_script) original_call <- simulation_results$call
  simulation_results <- simulation_results$results

  # Save the results per method
  method_names <- names(methods)
  for (i in 1:length(methods)) {

    tmp_results <- extract_named_elements(simulation_results, method_names[i])
    tmp_results <- dplyr::bind_rows(tmp_results)
    # results <- list("results" = tmp_results)
    results <- tmp_results

    if (mimic_script) {
      params <- sim_args
      method_info <- methods[[method_names[i]]]
      params[names(method_info)] <- method_info
      indexr::save_objects("./rds", results, args_list = params, overwrite = TRUE)
    } else {
      ## Clean up the call argument
      results$call <- original_call
      results$call["methods"] <- NULL
      method_info <- methods[[method_names[i]]]
      results$call[names(method_info)] <- method_info
      results$call[1] <- "sim" ## Need more elegant solution

      indexr::save_objects("./rds", results, overwrite = overwrite, get_script_name = FALSE)
    }

  }

}
#' Title
#'
#' @param methods
#' @param simulation_function
#' @param simulation_arguments
#' @param iterations
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
sim <- function(methods, simulation_function, simulation_arguments,
                iterations = 1000, seed = 1234) {

  ## Get original call
  original_call <- match.call()

  # Combine the method information and simulation information
  simulation_arguments <- list(
    "methods" = methods,
    "simulation_function" = simulation_function,
    "simulation_arguments" = simulation_arguments,
    "iterations" = iterations,
    "seed" = seed
  )

  # Run the simulation
  res <- lapply(1:iterations, function(x) do.call(sim_iteration, c(iteration = list(x), simulation_arguments)))

  return(list(call = original_call, results = res))

}
sim_parallel <- function(methods, simulation_function, simulation_arguments,
                         iterations = 1000, seed = 1234, nClust = 2,
                         libsEx, clusterEv, clusterEx) {

  ## Get original call, remove parallel arguments
  original_call <- match.call()
  original_call$nClust <- NULL
  original_call$libsEx <- NULL
  original_call$clusterEv <- NULL
  original_call$clusterEx <- NULL

  # Combine the method information and simulation information
  simulation_arguments <- list(
    "methods" = methods,
    "simulation_function" = simulation_function,
    "simulation_arguments" = simulation_arguments,
    "iterations" = iterations,
    "seed" = seed
  )

  ## Need to make more general for other users
  cl <- parallel::makeCluster(nClust)
  if (!is.null(clusterEx)) clusterExport(cl, clusterEx)
  if (!is.null(clusterEv)) clusterExport(cl, c("clusterEv", "libsEx"), envir = environment())
  clusterEvalQ(cl, {
    if(!is.null(clusterEv)) {clusterEv()}
    if(!is.null(libsEx)) lapply(libsEx, function(lib) {library(lib, character.only = TRUE)})
  })

  res <- parallel::parLapply(cl, 1:iterations, function(x) do.call(sim_iteration, c(iteration = list(x), simulation_arguments)))
  parallel::stopCluster(cl)

  return(list(call = original_call, results = res))

}
sim_iteration <- function(iteration, methods, simulation_function,
                          simulation_arguments = list(), iterations = 1000,
                          seed = 1234) {

  method_names <- names(methods)

  set.seed(seed)
  seeds <- round(runif(iterations) * 1e9)

  ret <- list()

  set.seed(seeds[iteration]) ## Could this make the methods correlated if ran separately?
  data <- do.call(simulation_function, simulation_arguments)

  for (j in 1:length(methods)) {

    time_taken <- system.time({
      results <- do.call(methods[[method_names[j]]]$method, c(list(X = data$X, y = data$y), methods[[method_names[j]]]$method_arguments))
    })
    ret[[method_names[j]]] <- results %>%
      mutate(truth = data$beta, iteration = iteration, method = method_names[j], time = time_taken["elapsed"])

  }

  return(ret)

}
# Function to extract elements associated with a specified name
extract_named_elements <- function(lst, name) {
  # Apply function to each sublist
  extracted <- lapply(lst, function(sublist) {
    # Check if the sublist is actually a list and contains the specified name
    if (is.list(sublist) && name %in% names(sublist)) {
      # Return the element associated with the specified name
      return(sublist[[name]])
    } else {
      return(NULL)
    }
  })

  # Filter out any NULL values (in case the name was not found in some sublists)
  extracted <- extracted[!sapply(extracted, is.null)]

  return(extracted)
}
