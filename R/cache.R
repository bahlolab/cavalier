
# environment to store cached tables, e.g. gtex expression, omim genemap etc
cavalier_cache <- new.env()

# execute a function and save to disk as filename unless filename exists, then load from file
# useful to cache downloaded files
#' @importFrom stringr str_ends str_c
cache <- function(fun, name,
                  disk = FALSE,
                  ver = '') 
{
  # attempt to return from memory
  value <- cavalier_cache[[name]]
  if (!is.null(value)) { return(value) }
  
  # attempt to return from disk
  if (disk) {
    cache_dir <- get_cavalier_opt('cache_dir')
    cache_file <- file.path(cache_dir, str_c(name, ver, '.rds'))
    if (file.exists(cache_file)) { 
      value <- readRDS(cache_file)
      cavalier_cache[[name]] <- value
      return(value)
    }
  }
  
  # evaluate function and save result
  value <- fun()
  cavalier_cache[[name]] <- value
  
  # save to disk
  if (disk) {
    if (!dir.exists(cache_dir)) { dir.create(cache_dir, recursive = TRUE) }
    # save to tempfile to avoid possible race condition
    tmp_fn <- tempfile(pattern = basename(cache_file) %>% str_c('.'),
                       tmpdir = cache_dir)
    saveRDS(value, tmp_fn)
    file.rename(tmp_fn, cache_file)
  }
  
  return(value)
}

#' @importFrom purrr walk map
#' @export
clear_cache <- function(mem = TRUE, disk = FALSE)
{
  if (disk) {
    cache_dir <- get_cavalier_opt('cache_dir')
    if (dir.exists(cache_dir)) {
      file.remove(list.files(cache_dir, full.names = TRUE))
    }
  }
  if (mem) {
    rm(list = ls(envir = cavalier_cache), envir = cavalier_cache)
  }
}
