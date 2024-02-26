rm(list = ls())
library(dplyr)
library(stringr)
library(uuid)
library(readr)

# Paths for old and new files
old_folder <- "/Users/loganharris/github/lasso-boot/rds"
new_folder <- "/Users/loganharris/github/lasso-boot/new_rds"
crosswalk_path <- file.path(new_folder, "crosswalk.csv")

# Read the existing crosswalk file or create a new one if it doesn't exist
if (file.exists(crosswalk_path)) {
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
} else {
  crosswalk <- data.frame(data = character(), n = integer(), rate = numeric(),
                          correlation_structure = character(), correlation = numeric(),
                          method = character(), ci_method = character(),
                          nominal_coverage = numeric(), p = integer(), uuid4 = character(),
                          timestamp = character(), stringsAsFactors = FALSE)
}

# List all relevant files in the old folder
datatype <- "abn"
file_list <- list.files(old_folder, pattern = glue("^{datatype}\\(\\d+\\)_"), full.names = TRUE)

for (file_path in file_list) {
  # Extract parameters from file name
  file_name <- str_remove(basename(file_path), ".rds")
  parts <- strsplit(file_name, "_")[[1]]
  modifier <- ifelse(length(parts) > 7, parts[8], NA)
  base_params <- list(data = datatype,
                      snr = as.numeric(str_extract(parts[2], "\\d+")),
                      # sd = as.numeric(gsub(glue("{datatype}\\((\\d+)\\).*"), "\\1", parts[1])),
                      # rate = as.numeric(gsub(glue("{datatype}\\((\\d+)\\).*"), "\\1", parts[1])),
                      a = as.numeric(gsub(glue("{datatype}\\((\\d+)\\).*"), "\\1", parts[1])),
                      b = 2,
                      correlation_structure = parts[3],
                      correlation = as.numeric(gsub("rho(.*)", "\\1", parts[4])),
                      method = parts[5],
                      ci_method = "quantile",
                      nominal_coverage = as.numeric(gsub("alpha(.*)", "\\1", parts[6])),
                      p = as.numeric(gsub("p(.*)", "\\1", parts[7])),
                      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))


  # Load the .RData file
  load(file_path)

  # Process each n value
  n_values <- unique(c(per_var$n, per_dataset$n))
  for (n in n_values) {
    # Generate a new UUID for each n
    uuid4 <- UUIDgenerate()

    # Create parameters list for this n
    params <- c(base_params, list(n = n, modifier = modifier, uuid4 = uuid4))

    # Filter data for this n
    per_var_n <- per_var %>% filter(n == n)
    per_dataset_n <- per_dataset %>% filter(n == n)

    # Save the filtered objects with the new UUID
    save_path <- file.path(new_folder, paste0(uuid4, ".rds"))
    save(per_var_n, per_dataset_n, file = save_path)

    # Create a new row for the crosswalk with specified column order
    column_order <- c("data", "n", "p", "snr", "rate", "a", "b", "sd",
                      "correlation_structure", "correlation", "correlation_noise",
                      "method", "ci_method", "nominal_coverage", "modifier", "uuid4", "timestamp")
    new_row <- setNames(nm = column_order, vector("list", length = length(column_order)))

    for (col in column_order) {
      new_row[[col]] <- ifelse(col %in% names(params), params[[col]], NA)
    }

    # Add entry to crosswalk
    crosswalk <- rbind(crosswalk, new_row)
  }
}

# Save the updated crosswalk
colnames(crosswalk)
write_csv(crosswalk, crosswalk_path, col_names = TRUE)
