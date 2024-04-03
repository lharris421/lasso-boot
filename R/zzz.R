.onLoad <- function(libname, pkgname) {
  # Set location for saving rds files
  res_dir <- switch(Sys.info()['user'],
                    'pbreheny' = '~/res/lasso-boot',
                    'loganharris' = '../lasso-boot')
  rds_path <<- glue::glue("{res_dir}/rds/")
  local_dir <- glue::glue("{res_dir}/local")

  unloadNamespace("hdrm")
  unloadNamespace("ncvreg")
  stopifnot(dir.exists(local_dir))
  .libPaths(local_dir)
  library(ncvreg)
}
