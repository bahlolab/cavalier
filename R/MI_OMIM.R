
get_mi_omim_version_latest <- function() {
  str_c(get_cavalier_opt("mi_omim_github_url"), 'releases') %>% 
    retry(verb = 'GET') %>% 
    content(encoding = 'UTF-8') %>% 
    rvest::html_nodes(".d-flex") %>%
    rvest::html_nodes(".d-inline")  %>%
    rvest::html_text(trim = TRUE) %>% 
    keep(str_detect, '^[0-9]{4}-[0-9]{2}-[0-9]{2}$') %>%
    sort() %>% 
    last()
}

get_mi_omim_version <- function(db_mode = get_cavalier_opt("database_mode"),
                                ver = get_cavalier_opt("mi_omim_version")) {
  if (is.null(ver)) {
    ver <-
      get_version(
        resource_name  = 'MI_OMIM',
        cache_name = 'mi_omim_names',
        cache_subdir = 'MI_OMIM',
        func_online = get_mi_omim_version_latest,
        db_mode = db_mode
      )
  }
  
  message("Using MI OMIM version ", ver)
  
  ver
}

get_mi_omim_names <- function(mi_omim_version = get_mi_omim_version()) {
  
  fun <- function() {
    mim_omim_url <- str_c(
      get_cavalier_opt("mi_omim_github_url"), 
      'releases/download/', 
      mi_omim_version,
      '/omim.sssom.tsv'
    )
    
    message('Downloading ', mim_omim_url)
    
    mim_omim_names <-
      retry('GET', mim_omim_url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(comment = '#', col_types = cols(.default = 'c')) %>% 
      select(omim_id = 1, omim_name = 2) %>% 
      distinct() %>% 
      arrange_all()
    
  }
  
  cache(
    fun = fun,
    name = 'mi_omim_names',
    subdir = 'MI_OMIM',
    ver = mi_omim_version,
  )
}


# combine HPO disease->gene mappings with MI_OMIM term names
#' @export
get_omim_disease_map <- function(gene_key = c('ensembl', 'entrez')) {
  
  get_genes_to_phenotype() %>% 
    select(entrez_id, symbol, disease_id, hpo_term_id) %>% 
    filter(str_starts(disease_id, 'OMIM:')) %>% 
    left_join(hpo_mendelian_inheritnace, by = 'hpo_term_id') %>% 
    group_by(entrez_id, symbol, disease_id) %>% 
    summarise(inheritance = str_c(sort(na.omit(inheritance)), collapse = '/'),
              .groups = 'drop') %>% 
    mutate(inheritance = if_else(nchar(inheritance) == 0, NA_character_, inheritance)) %>% 
    select(entrez_id, symbol, disease_id, inheritance) %>% 
    distinct() %>% 
    mutate(
      ensembl_id = coalesce(
        hgnc_entrez2ensembl(entrez_id),
        hgnc_sym2ensembl(symbol),
      ),
      .after = entrez_id
    ) %>% 
    left_join(
      get_mi_omim_names() %>% 
        select(disease_id = omim_id, disease_name = omim_name)
    ) %>% 
    mutate(disease_name = if_else(is.na(disease_name), '???', disease_name))
}

