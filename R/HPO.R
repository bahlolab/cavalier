
#' Mendelian inheritance terms for assigning inheritance to HPO diseases
hpo_mendelian_inheritnace <- tibble::tribble(
  ~ hpo_term_id, ~ inheritance,
  'HP:0000006' , 'AD',
  'HP:0000007' , 'AR',
  'HP:0001417' , 'XL',
  'HP:0001427' , 'MT'
)

get_hpo_version_latest <- function() {
  str_c(get_cavalier_opt("hpo_github_url"), 'tags') %>% 
    retry(verb = 'GET') %>% 
    content(encoding = 'UTF-8') %>% 
    rvest::html_nodes(".d-flex") %>%
    rvest::html_nodes(".d-inline")  %>%
    rvest::html_text(trim = TRUE) %>% 
    keep(str_detect, '^v[0-9]{4}-[0-9]{2}-[0-9]{2}$') %>%
    sort() %>% 
    last()
}

#' Get HPO version from GitHub or disk cache
get_hpo_version <- function(db_mode = get_cavalier_opt("database_mode"),
                            ver = get_cavalier_opt("hpo_version")
    ) 
{
  
  if (is.null(ver)) {
    ver <-
      get_version(
        resource_name  = 'HPO',
        cache_name = 'genes_to_phenotype.phenotype_to_genes',
        cache_subdir = 'HPO',
        func_online = get_hpo_version_latest,
        db_mode = db_mode
      )
  }
  
  message("Using HPO version ", ver)
  
  ver
}

#' Retrieve HPO genes_to_phenotype and phenotype_to_genes
#' 
#' Download or load from disk cache
get_hpo_g2p_p2g <- function(ver = get_hpo_version()) {
  
  
  fun <- function() {
    
    g2p_url <- str_c(
      get_cavalier_opt("hpo_github_url"), 
      'releases/download/', 
      ver,
      '/genes_to_phenotype.txt'
    )
    
    message('Downloading ', g2p_url)
    g2p <-
      retry('GET', g2p_url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('entrez_id',
                             'symbol',
                             'hpo_term_id',
                             'hpo_term_name', 
                             'frequency', 
                             'disease_id'),
               skip = 1,
               col_types = 'iccccc') %>% 
      select(-frequency)
    
    p2g_url <- str_c(
      get_cavalier_opt("hpo_github_url"),
      'releases/download/',
      ver,
      '/phenotype_to_genes.txt'
    )
    
    message('Downloading ', p2g_url)
    p2g <-
      retry('GET', p2g_url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('hpo_term_id',
                             'hpo_term_name',
                             'entrez_id', 
                             'symbol',
                             'disease_id'),
               col_types = 'ccicc',
               skip = 1)
    
    list(genes_to_phenotype = g2p,
         phenotype_to_genes = p2g)
    
  }
  
  cache(
    fun = fun,
    name = 'genes_to_phenotype.phenotype_to_genes',
    subdir = 'HPO',
    ver = ver,
  )
}

#' HPO genes_to_phenotype table
get_genes_to_phenotype <- function() 
{
  get_hpo_g2p_p2g()$genes_to_phenotype
}

#' HPO phenotype_to_genes table
get_phenotype_to_genes <- function() 
{
  get_hpo_g2p_p2g()$phenotype_to_genes
}


#' Simplified mapping from gene to phenotype from HPO
get_gene_disease_map <- function(source = c('ALL', 'OMIM', 'ORPHA')) 
{
  source <- match.arg(source)
  
  fun <- function() {
    get_genes_to_phenotype() %>% 
      select(entrez_id, symbol, disease_id, hpo_term_id) %>% 
      group_by(entrez_id, symbol, disease_id) %>% 
      left_join(hpo_mendelian_inheritnace, by = 'hpo_term_id') %>% 
      summarise(inheritance = str_c(sort(na.omit(inheritance)), collapse = '/'),
                .groups = 'drop') %>% 
      mutate(inheritance = if_else(nchar(inheritance) == 0, NA_character_, inheritance)) %>% 
      select(entrez_id, symbol, disease_id, inheritance) %>% 
      distinct()
  }
  
  result <-
    cache(
    fun = fun,
    name = 'gene_disease_map',
    subdir = 'HPO',
    ver = get_hpo_version()
  )
  
  if (source == 'ALL') {
    return(result)
  } else {
    return(
      filter(result, str_detect(disease_id, str_c('^', source)))
    )
  }
}

#' Mapping of hpo_term_id to hpo_term_name
get_hpo_term_names <- function(ver = get_hpo_version())
{
  bind_rows(
    get_phenotype_to_genes() %>%
      select(hpo_term_id, hpo_term_name),
    get_genes_to_phenotype() %>%
      select(hpo_term_id, hpo_term_name)
  ) %>% distinct()
}

#' convert hpo_term_id to hpo_term_name
hpo_id2name <- function(hpo_term_ids)
{
  term_names <- get_hpo_term_names()
  with(term_names, hpo_term_name[match(hpo_term_ids, hpo_term_id)])
}

#' Get a gene list from HPO phenotype_to_genes table
#' @export

get_hpo_gene_list <- function(hpo_id, prefer_omim = TRUE, hpo_version = get_hpo_version()) {
  
  assert_that(
    is_scalar_character(hpo_id),
    !is.na(hpo_id),
    str_detect(hpo_id, 'HP:\\d+$'),
    is_scalar_character(hpo_version) || is.null(hpo_version)
  )
  
  term_name <- hpo_id2name(hpo_id)
  
  get_phenotype_to_genes() %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_id, symbol, disease_id) %>% 
    distinct() %>% 
    left_join(
      get_gene_disease_map(source = 'ALL') %>%
        select(-symbol),
      by = c('entrez_id', 'disease_id')) %>% 
    group_by(entrez_id) %>% 
    filter(!prefer_omim | str_starts(disease_id, 'OMIM') | !any(str_starts(disease_id, 'OMIM'))) %>% 
    ungroup() %>% 
    mutate(symbol = coalesce(
      hgnc_entrez2sym(entrez_id),
      hgnc_sym2sym(symbol),
      symbol)) %>% 
    mutate(list_version = get_hpo_version()) %>% 
    mutate(list_id = hpo_id,
           list_name = term_name) %>% 
    select(list_id, list_name, list_version, symbol, entrez_id, disease_id, inheritance)
}
