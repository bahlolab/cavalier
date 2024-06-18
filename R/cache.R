
#' Save or load disk cache
#'
#' Store result of fun() in RDS file on disk, or load from disk if it exists
#' will store files in one of:
#'    <cache_dir>/<name>.rds
#'    <cache_dir>/<name>.ver_<version>.rds
#'    <cache_dir>/<subdir>/<name>.rds
#'    <cache_dir>/<subdir>/<name>.ver_<version>.rds
#'    
#' @importFrom stringr str_ends str_c
#' @importFrom dplyr last
cache <- function(name,
                  fun = NULL,
                  version = NULL,
                  subdir = NULL,
                  cache_dirs = get_cache_dirs())
{
  
  assert_that(
    is_scalar_character(name),
    is.function(fun) | is.null(fun),
    is_scalar_character(version) | is.null(version),
    is_scalar_character(subdir) | is.null(subdir),
    is.character(cache_dirs) & length(na.omit(cache_dirs)) > 0
  )
  
  name_ver <- `if`(!is.null(version), str_c(name, '.version_', version), name)
  
  # get cached object paths
  if (!is.null(subdir)) {
    paths <- file.path(cache_dirs, subdir, str_c(name_ver, '.rds'))
  } else {
    paths <- file.path(cache_dirs, str_c(name_ver, '.rds'))
  }
  
  # read first cached object if exists
  for (path in paths) {
    if (file.exists(path)) {
      message('Loading ', name_ver, ' from cache: ', path)
      return(readRDS(path))
    }
  }
  
  if (is.null(fun)) {
    stop("No cache exists for ", name_ver)
    return(NULL)
  }
  
  # Evaluate fun and attempt to store in first path
  value <- fun()
  for (path in paths) {
    dir <- dirname(path)
    tryCatch(
      {
        if (!dir.exists(dir)) { dir.create(dir, recursive = T) }
        tmp <- tempfile(pattern = str_c(name_ver, '.tmp_'), tmpdir = dir)
        saveRDS(value, file = tmp, compress = TRUE)
        file.rename(tmp, path)
        message("Saved ", name_ver, ' to ', path)
      },
      error = function(e) NULL
    )
    if (file.exists(path)) {
      return(value)
    }
  }

  warning('Could not save ', name_ver, ' to cache')

  return(value)
}

get_latest_cached_version <- function(
    name,
    subdir = NULL,
    cache_dirs = get_cache_dirs()) 
{
  
  if (!is.null(subdir)) {
    paths <- file.path(cache_dirs, subdir)
  } else {
    paths <- file.path(cache_dirs)
  }
  
  version <-
    list.files(paths) %>% 
    keep(str_detect, str_c('^', name, '\\.version_.+\\.rds$')) %>% 
    unique() %>% 
    str_extract('(?<=\\.version_).+(?=\\.rds$)') %>% 
    sort_versions() %>% 
    last()
  
  if (is.na(version)) {
    return(NULL)
  }
  
  return(version)
}

#' Get the latest version from online on cache depending on database_mode
get_version <- function(
    resource_name,
    cache_name,
    cache_subdir,
    func_online,
    db_mode
)
{
  version <- NULL
  
  version_cached <- get_latest_cached_version(
    name = cache_name,
    subdir = cache_subdir
  )
  
  if (db_mode == "offline" & !is.null(version_cached)) {
    message("Using ", resource_name, " version ", version_cached)
    return(version_cached)
  }
  
  version_latest <- tryCatch(
    func_online(),
    error = function(e) NULL
  )
  if (!is_scalar_character(version_latest) || is.na(version_latest)) {
    message("Failed to retrieve latest ", resource_name, " version")
    version_latest <- NULL
  }
  
  if (!is.null(version_latest)) {
    message("Using ", resource_name, " version ", version_latest)
    return(version_latest)
  } else if (db_mode == 'online') {
    stop("Failed to retrieve latest ",  resource_name, " version")
  }
  
  if (!is.null(version_cached)) {
    if (db_mode == 'fallback') {
      message('Falling back to cached ', resource_name, ' version')
    }
    message("Using ", resource_name, " version ", version_cached)
    return(version_cached)
  } else {
    stop("No cached version and failed to retrieve latest ",  resource_name, " version")
  }
  
}


#' Store or load extenal file to disk
cache_file <- function(url,
                       version = NULL,
                       subdir = NULL,
                       cache_dirs = get_cache_dirs())
{
  
  assert_that(
    is_scalar_character(version) | is.null(version),
    is_scalar_character(subdir) | is.null(subdir),
    is.character(cache_dirs) & length(na.omit(cache_dirs)) > 0
  )
  
  # get cached object paths
  if (!is.null(subdir)) {
    paths <- file.path(cache_dirs, subdir, basename(url))
  } else {
    paths <- file.path(cache_dirs, basename(url))
  }
  
  # read first cached object if exists
  for (path in paths) {
    if (file.exists(path)) {
      message('Using ', path)
      return(path)
    }
  }
  
  # Download file attempt to store in first path
  for (path in paths) {
    dir <- dirname(path)
    tryCatch(
      {
        if (!dir.exists(dir)) { dir.create(dir, recursive = T) }
        tmp <- tempfile(pattern = str_c(basename(url), '.tmp_'), tmpdir = dir)
        download.file(url, destfile = tmp)
        file.rename(tmp, path)
        message("Downloaded ", url, ' to ', path)
      },
      error = function(e) NULL
    )
    if (file.exists(path)) {
      return(path)
    }
  }
  
  stop('Could not download ', url)
}

get_cache_dirs <- function() {

  # system dir - useful for prebuilt cache in containers
  cache_dirs <- tryCatch(
    file.path(system.file(package = 'cavalier', mustWork = TRUE), 'cache'),
    error = function(e) character()
  )

  # user dir, defaults to '~/.cavalier_cache'
  user_cache_dir <- get_cavalier_opt('cache_dir')
  if (!is.null(user_cache_dir)) {
    # increment whenever breaking changes are made to cached files
    CACHE_VER <- "v1"
    user_cache_dir <- file.path(user_cache_dir, CACHE_VER)
    if (!file.exists(user_cache_dir)) {
      dir.create(user_cache_dir, recursive = TRUE)
    }
    cache_dirs <- c(user_cache_dir, cache_dirs)
  }

  assert_that(length(cache_dirs) > 0)

  return(cache_dirs)

}

#' Removes memory cache for functions memoised with memoise::memoise
#' 
#' @export
clear_memory_caches <- function() {
  getNamespace('cavalier') %>% 
    as.list() %>% 
    keep(memoise::is.memoised) %>% 
    walk(function(x) memoise::forget(x))
}

#' Build local disk caches
#' 
#' Runs all functions that use local disk caching to build a complete cache for running in offline or fallback mode
#' @export
build_caches <- function(
    GeVIR = TRUE,
    HGNC = TRUE,
    GTEx = TRUE,
    IGV = TRUE,
    UCSC = TRUE,
    PanelApp = TRUE,
    HPO = TRUE,
    HPO_disease_names = TRUE) 
{
  if (HGNC) {
    message('Building HGNC cache')
    invisible(get_hgnc_complete())
    message('HGNC done')
  }
  if (GeVIR) {
    message('Building GeVIR cache')
    invisible(get_gevir_table())
    message('GeVIR done')
  }
  if (GTEx) {
    message('Building GTEx cache')
    invisible(get_gtex_expression())
    message('GTEx done')
  }
  if (IGV) {
    message('Downloading IGV geneomes')
    invisible(get_igv_genome('hg38'))
    invisible(get_igv_genome('hg19'))
    message('IGV done')
  }
  if (UCSC) {
    message('Building UCSC assembly gaps')
    ref_genome <- get_cavalier_opt('ref_genome')
    set_cavalier_opt(ref_genome = 'hg38')
    invisible(get_centromeres_gaps())
    set_cavalier_opt(ref_genome = 'hg19')
    invisible(get_centromeres_gaps())
    set_cavalier_opt(ref_genome = ref_genome)
    message('UCSC assembly gaps done')
  }
  if (HPO) {
    message('Building HPO cache')
    invisible(get_genes_to_phenotype())
    invisible(get_phenotype_to_genes())
    invisible(get_gene_disease_map())
    message('HPO done')
  }
  if (HPO_disease_names) {
    message('Building HPO disease names cache')
    message('This may take some time...')
    invisible(build_disease_name_cache())
    message('HPO disease names cache done')
  }
  if (PanelApp) {
    message('Building PanelApp cache')
    message('This will take some time...')
    invisible(build_panelapp_cache(sources = names(get_cavalier_opt('panelapp_urls'))))
    message('PanelApp cache done')
  }
}

