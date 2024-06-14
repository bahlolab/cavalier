
# environment to store cached tables, e.g. gtex expression, omim genemap etc
cavalier_cache <- new.env()

# execute a function and save to disk as filename unless filename exists, then load from file
# useful to cache downloaded files
#' @importFrom stringr str_ends str_c
cache <- function(fun,
                  name,
                  disk = FALSE,
                  overwrite = FALSE,
                  ver = NULL,
                  subdir = NULL,
                  store = TRUE) 
{
  # attempt to return from memory
  value <- cavalier_cache[[name]]
  if (!is.null(value) && !overwrite) { return(value) }
  
  # attempt to return from disk
  if (disk) {
    cache_dir <- get_cache_dir(subdir = subdir)
    cache_file <- `if`(!is.null(ver),
                       file.path(cache_dir, str_c(name, '.ver_', ver, '.rds')),
                       file.path(cache_dir, str_c(name, '.rds')))
    if (file.exists(cache_file) & !overwrite) { 
      value <- readRDS(cache_file)
      cavalier_cache[[name]] <- value
      return(value)
    }
  }
  
  # evaluate function and save result
  value <- fun()
  
  if (store) {
    cavalier_cache[[name]] <- value
    
    # save to disk
    if (disk) {
      if (!dir.exists(cache_dir)) { dir.create(cache_dir, recursive = TRUE) }
      # save to tempfile to avoid possible race condition
      tmp_fn <- tempfile(pattern = basename(cache_file) %>% str_c('.'),
                         tmpdir = cache_dir)
      saveRDS(value, tmp_fn)
      # attempt to avoid race conditions
      if (file.exists(cache_file) && !overwrite) {
        file.remove(tmp_fn)
      } else {
        file.rename(tmp_fn, cache_file)
      }
      
    }
  }
  
  return(value)
}

#' @importFrom purrr walk map
#' @export
clear_cache <- function(mem = TRUE, disk = FALSE)
{
  if (disk) {
    cache_dir <- get_cache_dir()
    if (dir.exists(cache_dir)) {
      file.remove(list.files(cache_dir, full.names = TRUE))
    }
  }
  if (mem) {
    rm(list = ls(envir = cavalier_cache), envir = cavalier_cache)
  }
}

get_cache_dir <- function(subdir = NULL) {
  
  # update whenever breaking changes are made to cached files
  CACHE_VER <- "0.2" 
  
  dir <- file.path(get_cavalier_opt('cache_dir'), CACHE_VER )
  
  if (!is.null(subdir)) {
    dir <- file.path(dir, subdir)
  }
  
  return(dir)
}

build_caches <- function(
    GeVIR = TRUE,
    HGNC = TRUE,
    GTEx = TRUE,
    HPO = TRUE,
    PanelApp = TRUE,
    OMIM_disease_names = TRUE) 
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
  if (HPO) {
    message('Building HPO cache')
    invisible(get_genes_to_phenotype())
    invisible(get_phenotype_to_genes())
    invisible(get_gene_disease_map(source = 'OMIM'))
    invisible(get_gene_disease_map(source = 'ORPHA'))
    message('HPO done')
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
  if (PanelApp) {
    message('Building PanelApp cache')
    message('This will take some time...')
    invisible(build_panelapp_cache(sources = c('PAA', 'PAE')))
    message('PanelApp disease name done')
  }
  if (OMIM_disease_names) {
    message('Building OMIM disease names cache')
    message('This will take some time...')
    invisible(build_disease_name_cache())
    message('OMIM disease name done')
  }
}
