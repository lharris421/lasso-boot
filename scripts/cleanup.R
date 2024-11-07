source("./scripts/setup/setup.R")
indexr::start_tagging(rds_path)
## Run manuscript
indexr::cleanup(rds_path)
indexr::close_tagging(rds_path)
