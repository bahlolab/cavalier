
.onLoad <- function(libname, pkgname) {
  
  # create memory-cached version of cavalier functions, most useful for api requests

  hpo_api_get <<- memoise::memoise(hpo_api_get)
  
  
  invisible()
}