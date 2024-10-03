library(digest)

# Read the crosswalk CSV
crosswalk_path <- glue("{rds_path}/crosswalk.csv")
crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE, na.strings = "")

# Create a vector to store hashes
hashes <- character(nrow(crosswalk))

# Iterate over each row of the crosswalk
for (i in 1:nrow(crosswalk)) {
  # Extract the row, ignore last two columns (uuid4 and timestamp)
  row <- crosswalk[i, 1:(ncol(crosswalk) - 2)]

  # Convert the row to a list and retain names
  args_list <- as.list(row)

  # Remove NA values while preserving names
  args_list <- Filter(function(x) !is.na(x), args_list)

  # Generate hash and store it
  hash <- digest(args_list, algo = "xxhash64")
  hashes[i] <- hash
}

# Find duplicate hashes
dup_hashes <- hashes[duplicated(hashes)]
if (length(dup_hashes) > 0) {
  # Identify rows causing duplicate hashes
  for (hash in dup_hashes) {
    dup_rows <- which(hashes == hash)
    cat("Duplicate hash:", hash, "\nRows:", toString(dup_rows), "\n\n")
  }
} else {
  cat("No duplicate hashes found.")
}
