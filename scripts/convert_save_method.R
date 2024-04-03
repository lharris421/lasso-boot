library(digest)
library(glue)

# Read the crosswalk CSV
crosswalk_path <- glue("../lasso-boot/archive_rds/crosswalk.csv")
crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# Function to determine and convert the data type of a column
convert_column <- function(column) {
  if (all(grepl("^[0-9]+$", column))) {
    return(as.numeric(column))
  } else {
    return(as.character(column))
  }
}

# Iterate over each row of the crosswalk
for (i in 1:nrow(crosswalk)) {
  # Extract the row, keep timestamp separate
  row <- crosswalk[i, 1:(ncol(crosswalk) - 2)]
  timestamp <- crosswalk$timestamp[i]

  # Apply the data type conversion function to each column
  row <- lapply(row, convert_column)

  # Convert the row to a named list and sort alphabetically by names
  args_list <- setNames(as.list(row), names(row))
  args_list <- args_list[order(names(args_list))]

  # Remove NA values while preserving names
  args_list <- Filter(function(x) !is.na(x), args_list)

  # Exclude timestamp from hash calculation
  args_list_no_timestamp <- args_list
  args_list_no_timestamp$timestamp <- NULL

  # Generate new filename using hash
  hash <- digest(args_list_no_timestamp, algo = "xxhash64")
  new_file_path <- glue("../lasso-boot/rds/{hash}.rds")

  # Add the timestamp to args_list
  args_list$timestamp <- timestamp

  # Read the old RDS file
  old_uuid <- crosswalk$uuid4[i]
  old_file_path <- glue("../lasso-boot/archive_rds/{old_uuid}.rds")

  if (file.exists(old_file_path)) {
    # Create a new environment to load the objects
    new_env <- new.env()

    # Load old RDS file into the new environment
    load(old_file_path, envir = new_env)

    # Get the names of all objects loaded
    object_names <- ls(envir = new_env)

    # Prepare the list of objects to be saved, including args_list
    objects_to_save <- c(list(args_list = args_list), mget(object_names, envir = new_env))

    # Save all objects in the new file
    save(list = names(objects_to_save), file = new_file_path, envir = new_env)
  } else {
    warning(paste0("Old file not found for UUID: ", old_uuid))
  }

}
