
#' Get HPO version from GitHub or disk cache
#' @importFrom purrr keep
get_g4e_version <- function(db_mode = get_cavalier_opt("database_mode")) {
  
  func_online <- function() {
    str_c(get_cavalier_opt("g4e_github_url"), 'releases') %>% 
      retry(verb = 'GET') %>% 
      content(encoding = 'UTF-8') %>% 
      rvest::html_nodes(".d-flex") %>%
      rvest::html_nodes(".d-inline")  %>%
      rvest::html_text(trim = TRUE) %>% 
      keep(str_detect, '_v20[0-9]{2}-[0-9]{2}$') %>% 
      str_extract('v20[0-9]{2}-[0-9]{2}') %>% 
      sort() %>%
      last()
  }
  
  get_version(
    resource_name  = 'Genes4Epilepsy',
    cache_name = 'epilepsy_genes',
    cache_subdir = 'Genes4Epilepsy',
    func_online = func_online,
    db_mode = db_mode
  )
}

#' Retrieve Genes4Epilepsy from GitHub
#' 
#' https://github.com/bahlolab/Genes4Epilepsy
#' @importFrom dplyr across
get_g4e_full_list <- function(version = get_g4e_version()) {
  
  fun <- function() {
    url <-
      get_cavalier_opt("g4e_github_url") %>% 
      str_replace('github', 'raw.githubusercontent') %>% 
      str_c('main/EpilepsyGenes_', version, '.tsv')
    
    read_tsv(url) %>% 
      rename(
        hgnc_id = HGNC_ID,
        symbol = Gene,
        ensembl_gene_id = Ensemble_ID,
        entrez_id = Entrez_ID,
        omim_id = OMIM_ID ,
        inheritance = Inheritance,
        phenotype = `Phenotype(s)`
      ) %>% 
      mutate(across(where(is.numeric), as.integer))
  }
 
  cache(
    fun = fun,
    name = 'epilepsy_genes',
    version = version,
    subdir = 'Genes4Epilepsy'
  )
  
}

#' Get sublist of Genes4Epilepsy by phenotype
get_g4e_phenotype_list <- function(
    version = get_g4e_version(),
    phenotype = c('ALL', 'DEE', 'Focal', 'MCD', 'GGE', 'PME')
) 
{
  phenotype <- match.arg(phenotype)
  
  if (is.null(version)) {
    version <- get_g4e_version()
  }
  
  get_g4e_full_list(version) %>% 
    mutate(phenotype = str_c('ALL, ', phenotype)) %>% 
    separate_rows(phenotype) %>% 
    filter(phenotype == !!phenotype) %>% 
    select(-phenotype) %>% 
    mutate(
      list_id = str_c('G4E:', phenotype),
      list_version = version,
      .before = 1)
  
}