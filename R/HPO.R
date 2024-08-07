
#' Mendelian inheritance terms for assigning inheritance to HPO diseases
hpo_mendelian_inheritnace <- tibble::tribble(
  ~ hpo_term_id, ~ inheritance,
  'HP:0000006' , 'AD',
  'HP:0000007' , 'AR',
  'HP:0001417' , 'XL',
  'HP:0001427' , 'MT'
)

#' Get HPO version from GitHub or disk cache
#' @importFrom purrr keep
get_hpo_version <- function(db_mode = get_cavalier_opt("database_mode")) {
  
  func_online <- function() {
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
  
  get_version(
    resource_name  = 'HPO',
    cache_name = 'genes_to_phenotype.phenotype_to_genes',
    cache_subdir = 'HPO',
    func_online = func_online,
    db_mode = db_mode
  )
} 

#' Retrieve HPO genes_to_phenotype and phenotype_to_genes
#' 
#' Download or load from disk cache
get_hpo_g2p_p2g <- function(hpo_version = get_hpo_version()) {
  
  fun <- function() {
    
    g2p_url <- str_c(
      get_cavalier_opt("hpo_github_url"), 
      'releases/download/', 
      hpo_version,
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
      hpo_version,
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
    ver = hpo_version,
  )
}

#' HPO genes_to_phenotype table
get_genes_to_phenotype <- function(hpo_version = get_hpo_version()) 
{
  get_hpo_g2p_p2g(hpo_version)$genes_to_phenotype
}

#' HPO phenotype_to_genes table
get_phenotype_to_genes <- function(hpo_version = get_hpo_version()) 
{
  get_hpo_g2p_p2g(hpo_version)$phenotype_to_genes
}

#' Simplified mapping from gene to phenotype from HPO
get_gene_disease_map <- function(source = c('ALL', 'OMIM', 'ORPHA'),
                                 hpo_version = get_hpo_version()) 
{
  source <- match.arg(source)
  
  fun <- function() {
    get_genes_to_phenotype(hpo_version) %>% 
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
    ver = hpo_version
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
get_hpo_term_names <- function(hpo_version = get_hpo_version())
{
  bind_rows(
    get_phenotype_to_genes(hpo_version) %>%
      select(hpo_term_id, hpo_term_name),
    get_genes_to_phenotype(hpo_version) %>%
      select(hpo_term_id, hpo_term_name)
  ) %>% distinct()
}

#' convert hpo_term_id to hpo_term_name
hpo_id2name <- function(hpo_term_ids, hpo_version = get_hpo_version())
{
  term_names <- get_hpo_term_names(hpo_version)
  with(term_names, hpo_term_name[match(hpo_term_ids, hpo_term_id)])
}

#' Get a gene list from HPO phenotype_to_genes table
#' @export
#' @importFrom dplyr inner_join anti_join add_row
get_hpo_gene_list <- function(hpo_id, prefer_omim = TRUE, hpo_version = get_hpo_version()) {
  
  assert_that(
    is_scalar_character(hpo_id),
    !is.na(hpo_id),
    str_detect(hpo_id, 'HP:\\d+$'),
    is_scalar_character(hpo_version) || is.null(hpo_version)
  )
  
  if (is.null(hpo_version)) {
    hpo_version <- get_hpo_version()
  }
  
  term_name <- hpo_id2name(hpo_id, hpo_version)
  
  get_phenotype_to_genes(hpo_version) %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_id, symbol, disease_id) %>% 
    distinct() %>% 
    left_join(
      get_gene_disease_map(source = 'ALL', hpo_version = hpo_version) %>%
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

#' Query HPO API
#'@importFrom httr GET accept_json content
#'@importFrom stringr str_c str_replace_all str_remove
hpo_api_get <- function(query, prefix = '', suffix = '')
{
  
  url <- str_c(
    get_cavalier_opt("hpo_api_base_url"),
    prefix,
    query, 
    suffix,
    sep = '/') %>% 
    str_replace_all('(?<!:)\\/+', '/') %>% 
    str_remove('\\/$')
  
  response <- retry('GET', url, accept_json())
  
  if (response$status_code != 200) {
    stop(url, ' returned code ', response$status_code)
  }
  
  return(content(response))
}

#' Query HPO API for disease name with disease id
hpo_api_get_disease_names <- function(disease_ids) {
  
  assert_that(
    is.character(disease_ids),
    all(na.omit(str_detect(disease_ids, '^(OMIM|ORPHA):')))
  )
  
  num_failuers <- 0L
  
  disease_names <-
    map_chr(disease_ids, function(disease_id) {
      if (!is.na(disease_id) & num_failuers < get_cavalier_opt('hpo_api_max_failuers')) {
        result <- tryCatch(
          hpo_api_get(disease_id, prefix = 'network/annotation'),
          error = function(e) { 'ERROR' }
        )
        if (is_scalar_character(result) && result == 'ERROR') {
          num_failuers <- num_failuers + 1L
        }
        tryCatch(result$disease$name, error = function(e) NA_character_)
      } else {
        NA_character_
      } 
    })
  
  
  if (num_failuers > get_cavalier_opt('hpo_api_max_failuers')) {
    warning('exceeded ', get_cavalier_opt('hpo_api_max_failuers'), ' HPO API failures')
  }
  
  disease_names
}

#' Query HPO API for diseases associates with gene, returning table of ids and name
hpo_api_get_diseases_by_gene <- function(entrez_ids) {
  
  assert_that(is_integer(entrez_ids))
  
  num_failuers <- 0L
  
  disease_tbl <-
    map_df(entrez_ids, function(entrez_id) {
      value <-
        tibble(entrez_id,
               disease_id = NA_character_, 
               disease_name = NA_character_,
               mondo_id = NA_character_)
      
      if (!is.na(entrez_id) & num_failuers < get_cavalier_opt('hpo_api_max_failuers')) {
        result <- tryCatch(
          hpo_api_get(str_c('NCBIGene:', entrez_id), prefix = 'network/annotation'),
          error = function(e) { 'ERROR' }
        )
        if (is_scalar_character(result) && result == 'ERROR') {
          num_failuers <- num_failuers + 1L
        }
        value <- 
          tryCatch(
            bind_rows(result$diseases) %>% 
              mutate(entrez_id = entrez_id) %>% 
              select(entrez_id, 
                     disease_id = id, 
                     disease_name = name,
                     mondo_id = mondoId), 
            error = function(e) value)
      }
      
      return(value)
      
    }) %>% 
    mutate(across(where(is.character), ~ replace(., nchar(.) == 0, NA_character_)))
  
  if (num_failuers > get_cavalier_opt('hpo_api_max_failuers')) {
    warning('exceeded ', get_cavalier_opt('hpo_api_max_failuers'), ' HPO API failures')
  }
  
  return(disease_tbl)

}

#' Build a complete cache of disease names
#' 
#' Use as a fallback method if hpo_api_get_disease_names fails
#' @export
build_disease_name_cache <- function(version = get_hpo_version()) 
{

  fun <- function() {
    disease_name_tbl <- 
      get_gene_disease_map(source = 'ALL') %>% 
      select(entrez_id, disease_id) %>% 
      distinct() %>% 
      add_count(entrez_id, name = 'n_disease') %>% 
      group_by(disease_id) %>% 
      slice(which.max(n_disease)) %>% 
      ungroup() %>% 
      pull(entrez_id) %>% 
      unique() %>% 
      hpo_api_get_diseases_by_gene() %>% 
      select(disease_id, disease_name) %>% 
      distinct() %>% 
      arrange(disease_id)
  }
  
  result <- cache(
    fun = fun,
    name = 'disease_names',
    version = version,
    subdir = 'HPO'
  )
  
  if (memoise::is.memoised(get_disease_name_cache)) { 
    memoise::drop_cache(get_disease_name_cache) 
  }
  
  invisible(result)
}

#' Get latest disease name cache from disk
get_disease_name_cache <- function()
{
  ver <- get_latest_cached_version(name = 'disease_names', subdir = 'HPO')
  if (is_scalar_character(ver) && !is.na(ver)) {
    cache(
      name = 'disease_names', 
      subdir = 'HPO', 
      ver = ver
    )
  } else {
    tibble(disease_id = character(), disease_name = character())
  }
}


#' Convert HPO disease_id to disease name
#' 
#' Uses HPO API or disk cache depending on cavalier_opts$database_mode
get_disease_names <- function(
    disease_ids,
    db_mode = get_cavalier_opt("database_mode")
) 
{
  
  disease_name_tbl <- 
    tibble(disease_id = disease_ids,
           disease_name = NA_character_)
  
  # get disease names from API
  if (db_mode %in% c("fallback", "online")) {
    disease_name_tbl <- 
      disease_name_tbl %>% 
      mutate(disease_name = hpo_api_get_disease_names(disease_id))
  }
  
  # replace missing disease names with latest cached values
  if (db_mode %in% c("offline", "fallback")) {
    disease_name_tbl <-
      disease_name_tbl %>%
      left_join(get_disease_name_cache(),
                by = 'disease_id',
                suffix = c('', '.cache')) %>% 
      mutate(disease_name = coalesce(disease_name, disease_name.cache)) %>% 
      select(-disease_name.cache)
  }
  
  disease_name_tbl %>% 
    pull(disease_name)
}




