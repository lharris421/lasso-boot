library(uuid)
library(dplyr)
library(tidyr)

save_objects <- function(folder, args_list, overwrite = FALSE, ...) {
  # Read the crosswalk CSV
  crosswalk_path <- file.path(folder, "crosswalk.csv")
  ## May want to consider intialization later
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

  # Validate arguments
  valid_cols <- colnames(crosswalk)
  invalid_args <- setdiff(names(args_list), valid_cols)
  if (length(invalid_args) > 0) {
    stop("Invalid arguments provided: ", paste(invalid_args, collapse = ", "))
  }

  # Check for existing parameters
  existing <- check_parameters_existence(folder, args_list, check_for = "existing")
  if (!is.null(existing)) {
    if (overwrite) {
      message("Overwriting existing record.")
      uuid4 <- existing$uuid4[1] # Assuming the first match is the one to overwrite
    } else {
      warning("Existing record found. Set 'overwrite = TRUE' to overwrite.")
      return(invisible())
    }
  } else {
    # Generate UUIDv4
    uuid4 <- UUIDgenerate()
  }

  # Generate UUIDv4 and timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  args_list$uuid4 <- uuid4
  args_list$timestamp <- timestamp

  # Update or append the new row to the crosswalk
  if (overwrite) {
    crosswalk <- crosswalk[crosswalk$uuid4 != uuid4, ] # Remove old record
  }

  # Complete the missing columns with NA
  missing_cols <- setdiff(valid_cols, names(args_list))
  args_list[missing_cols] <- NA

  # Append the new row to the crosswalk and save
  crosswalk <- rbind(crosswalk, args_list)
  write_csv(crosswalk, crosswalk_path, col_names = TRUE)

  # Save objects to RDS
  save(..., file = file.path(folder, paste0(uuid4, ".rds")))
}

read_objects <- function(folder, params_grid) {

  check_parameters_existence(folder, params_grid, check_for = "missing", halt = TRUE)

  # Read the crosswalk CSV
  crosswalk_path <- file.path(folder, "crosswalk.csv")
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

  # Identify columns in crosswalk that are not all NA and are in params_list
  valid_cols <- names(crosswalk)[colSums(!is.na(crosswalk)) > 0 & names(crosswalk) %in% names(params_grid)]
  # Expand the parameters to all combinations
  params_combinations <- params_grid[valid_cols]

  # Perform an inner join to find matching rows in crosswalk
  uuids <- inner_join(crosswalk, params_combinations, by = names(params_combinations)) %>%
    pull(uuid4)
  print(uuids)

  # Read the corresponding .rds files and bind them
  for (i in 1:length(uuids)) {
    file_path <- file.path(folder, paste0(uuids[i], ".rds"))
    load(file_path, envir = globalenv())
  }
}

check_parameters_existence <- function(folder, params_list, check_for = c("missing", "existing"), halt = FALSE) {
  # Validate 'check_for' option
  check_for <- match.arg(check_for)

  # Read the crosswalk CSV
  crosswalk_path <- file.path(folder, "crosswalk.csv")
  if (!file.exists(crosswalk_path)) {
    stop("Crosswalk file does not exist at the specified path: ", crosswalk_path)
  }
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

  # Expand the parameters to all combinations (if they contain vectors)
  params_combinations <- expand.grid(params_list)

  # Convert params_combinations to a data frame for joining
  params_combinations <- as.data.frame(params_combinations, stringsAsFactors = FALSE)

  # Check for missing or existing parameters
  result <- switch(check_for,
                   missing = {
                     # Perform an anti join to find missing combinations
                     missing_combinations <- anti_join(params_combinations, crosswalk,
                                                       by = intersect(names(params_combinations), names(crosswalk)))
                     if (nrow(missing_combinations) > 0) {
                       if (halt) {
                         print(missing_combinations)
                         stop("Missing parameters found.")
                       }
                       missing_combinations
                     } else NULL
                   },
                   existing = {
                     # Perform an inner join to find existing combinations
                     existing_combinations <- inner_join(params_combinations, crosswalk,
                                                         by = intersect(names(params_combinations), names(crosswalk)))
                     if (nrow(existing_combinations) > 0) {
                       if (halt) {
                         print(existing_combinations)
                         stop("Existing parameters found.")
                       }
                       existing_combinations
                     } else NULL
                   }
  )

  return(result)
}
