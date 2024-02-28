library(stringr)


# Paths for old and new files
old_folder <- "/Users/loganharris/github/lasso-boot/rds"
new_folder <- "/Users/loganharris/github/lasso-boot/new_rds"
crosswalk_path <- file.path(new_folder, "crosswalk.csv")

# Read the existing crosswalk file or create a new one if it doesn't exist
if (file.exists(crosswalk_path)) {
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
}

# List all relevant files in the old folder
datatype <-"across"
file_list <- list.files(old_folder, pattern = glue("^{datatype}_"), full.names = TRUE)

for (file_path in file_list) {
  # Extract parameters from file name
  load(file_path)
  # load(str_remove(file_path, "_example"))
  file_name <- str_remove(basename(file_path), ".rds")
  parts <- strsplit(file_name, "_")[[1]]
  base_params <- list(data = "laplace",
                      rate = 2,
                      snr = 1,
                      n = 100,
                      p = 100,
                      correlation_structure = "exchangeable",
                      correlation = 0,
                      correlation_noise = NA,
                      method = "zerosample2",
                      ci_method = "quantile",
                      nominal_coverage = .2 * 100)

  save_objects(folder = new_folder, args_list = base_params, res, lambdas)
}
