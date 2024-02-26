library(uuid)
library(dplyr)
library(tidyr)

save_objects <- function(folder, args_list, ...) {
  # Read the crosswalk CSV
  crosswalk_path <- file.path(folder, "crosswalk.csv")
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

  # Validate arguments
  valid_cols <- colnames(crosswalk)
  invalid_args <- setdiff(names(args_list), valid_cols)
  if (length(invalid_args) > 0) {
    stop("Invalid arguments provided: ", paste(invalid_args, collapse = ", "))
  }

  # Generate UUIDv4 and timestamp
  uuid4 <- UUIDgenerate()
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  # Add UUIDv4 and timestamp to the args_list
  args_list$uuid4 <- uuid4
  args_list$timestamp <- timestamp

  # Complete the missing columns with NA
  missing_cols <- setdiff(valid_cols, names(args_list))
  args_list[missing_cols] <- NA

  # Append the new row to the crosswalk and save
  crosswalk <- rbind(crosswalk, args_list)
  write_csv(crosswalk, crosswalk_path, col_names = TRUE)

  # Save objects to RDS
  save(..., file = file.path(folder, paste0(uuid4, ".rds")))
}

read_objects <- function(folder, params_list) {
  # Read the crosswalk CSV
  crosswalk_path <- file.path(folder, "crosswalk.csv")
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

  # Expand the parameters to all combinations
  params_combinations <- expand.grid(params_list)

  # Filter the crosswalk based on the combinations
  filtered_crosswalk <- crosswalk
  for (param_name in names(params_combinations)) {
    param_values <- unique(params_combinations[[param_name]])
    if (param_name %in% names(crosswalk)) {
      filtered_crosswalk <- filtered_crosswalk[filtered_crosswalk[[param_name]] %in% param_values, ]
    } else {
      warning("Parameter not found in crosswalk: ", param_name)
    }
  }

  # Read the corresponding .rds files and bind them
  all_data <- lapply(filtered_crosswalk$uuid4, function(uuid) {
    file_path <- file.path(folder, paste0(uuid, ".rds"))
    if (file.exists(file_path)) {
      load(file_path)  # Changed from readRDS to load
      list(per_var_all, per_dataset_all)  # Assuming these are the common objects in your RDS files
    } else {
      warning("File not found: ", file_path)
      NULL
    }
  })

  # Optionally, bind the data together if needed
  # This step depends on the structure of your .rds files and how you want to combine them
  # Example for row-binding if your objects are data frames
  all_data_combined <- do.call(rbind, lapply(all_data, function(x) x$per_var_all))
  # Repeat for per_dataset_all or any other objects as needed

  return(all_data_combined)
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
