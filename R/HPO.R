
hpo_mendelian_inheritnace <- tibble::tribble(
  ~ hpo_term_id, ~ inheritance,
  'HP:0000006' , 'AD',
  'HP:0000007' , 'AR',
  'HP:0001417' , 'XL',
  'HP:0001427' , 'MT'
)

#' @importFrom purrr keep
latest_hpo_release <- function() {
  
  fun <- function() {
    # attempt to get latest hpo release from github, but otherwise use cached version in case server is down
    url <- get_cavalier_opt("hpo_github_tag_url")
    
    tryCatch(
      content(retry(url, verb = 'GET')) %>% 
        map_chr('name') %>% 
        keep(str_detect, 'v[0-9]{4}-[0-9]{2}-[0-9]{2}') %>% 
        sort(decreasing = TRUE) %>% 
        first() %>% 
        (function(x) { assert_that(is_scalar_character(x) & !is.na(x)); x }),
      error = function(e) {
        hpo_files <- list.files(get_cache_dir(), pattern = '^hpo\\..*\\v[0-9]{4}-[0-9]{2}-[0-9]{2}\\.rds$')
        if (length(hpo_files)) {
          ver <- 
            str_extract(hpo_files, 'v[0-9]{4}-[0-9]{2}-[0-9]{2}(?=\\.rds)') %>% 
            sort(decreasing = TRUE) %>% 
            first()
          warning("Coudn't retrieve latest HPO release from: ", tag_url, '. ',
                  "Using cached version ", ver, '.')
          ver
        } else {
          stop('could not get latest hpo release')
        }
      })
  } 
  
  cache(fun = fun, name = 'latest_hpo_release')
} 

# this contains childmost term for a phenotype heirachy
# use for annotating genes
get_genes_to_phenotype <- function() 
{

  fun <- function() {
    
    url <- str_c(get_cavalier_opt("hpo_github_base_url"), 'releases/download/', hpo_release, '/genes_to_phenotype.txt')
    
    retry('GET', url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('entrez_id',
                             'symbol',
                             'hpo_term_id',
                             'hpo_term_name', 
                             'frequency', 
                             'disease_id'),
               skip = 1,
               col_types = 'iccccc')
  }
  
  cache(
    fun = fun,
    name = 'genes_to_phenotype',
    subdir = 'HPO',
    ver = latest_hpo_release(),
    disk = TRUE
  )
}

# this contains all hpo terms
# use for creating gene lists from hpo terms
get_phenotype_to_genes <- function() 
{

  fun <- function() {
    
    url <- str_c(get_cavalier_opt("hpo_github_base_url"), 'releases/download/', hpo_release, '/phenotype_to_genes.txt')
    
    retry('GET', url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('hpo_term_id',
                             'hpo_term_name',
                             'entrez_id', 
                             'symbol',
                             'disease_id'),
               col_types = 'ccicc',
               skip = 1)
  }
  
  cache(
    fun = fun,
    name = 'phenotype_to_genes',
    subdir = 'HPO',
    ver = latest_hpo_release(),
    disk = TRUE
  )
}

get_gene_disease_map <- function(source = c('OMIM', 'ORPHA')) 
{
  source <- match.arg(source)
  
  fun <- function() {
    get_genes_to_phenotype() %>% 
      select(entrez_id, symbol, disease_id, hpo_term_id) %>% 
      filter(str_starts(disease_id, str_c(source, ':'))) %>% 
      group_by(entrez_id, symbol, disease_id) %>% 
      left_join(hpo_mendelian_inheritnace, by = 'hpo_term_id') %>% 
      summarise(inheritance = str_c(sort(na.omit(inheritance)), collapse = '/'),
                .groups = 'drop') %>% 
      select(entrez_id, symbol, disease_id, inheritance) %>% 
      distinct()
  }
  
  cache(
    fun = fun,
    name = str_c('gene_disease_map_', source),
    subdir = 'HPO',
    disk = TRUE,
    ver = latest_hpo_release()
  )
}

hpo_term_names <- function(hpo_term_ids)
{
  fun <- function() {
    distinct(
      bind_rows(
      get_phenotype_to_genes() %>%
        select(hpo_term_id, hpo_term_name),
      get_genes_to_phenotype() %>%
        select(hpo_term_id, hpo_term_name)
    ))
  }
  
  term_names <- cache(fun = fun, name = 'hpo_term_names')
  
  with(term_names, hpo_term_name[match(hpo_term_ids, hpo_term_id)])
}

#' @export
#' @importFrom dplyr inner_join anti_join add_row
get_hpo_gene_list <- function(hpo_id, prefer_omim = TRUE) {
  
  assert_that(
    is_scalar_character(hpo_id),
    !is.na(hpo_id),
    str_detect(hpo_id, 'HP:\\d+$')
  )
  
  term_name <- hpo_term_names(hpo_id)
  
  get_phenotype_to_genes() %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_id, symbol, disease_id) %>% 
    distinct() %>% 
    left_join(
      bind_rows(
        get_gene_disease_map(source = 'OMIM'),
        get_gene_disease_map(source = 'ORPHA')
      ) %>% select(-symbol),
      by = c('entrez_id', 'disease_id')) %>% 
    group_by(entrez_id) %>% 
    filter(!prefer_omim | str_starts(disease_id, 'OMIM') | !any(str_starts(disease_id, 'OMIM'))) %>% 
    ungroup() %>% 
    mutate(gene = coalesce(
      hgnc_entrez2sym(entrez_id),
      hgnc_sym2sym(symbol),
      symbol)) %>% 
    select(gene, disease_id,  inheritance) %>%
    arrange(gene, disease_id) %>% 
    mutate(version = latest_hpo_release()) %>% 
    mutate(list_id = hpo_id,
           list_name = term_name) %>% 
    select(list_id, list_name, everything())
}

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
  
  if (response$status_code == 200) {
    return(content(response))
  }
  return(NULL)
}

hpo_api_get_disease_names <- function(disease_ids) {
  
  assert_that(
    is.character(disease_ids),
    all(str_detect(disease_ids, '^(OMIM|ORPHA):'))
  )
  
  num_failuers <- 0L
  
  disease_names <-
    map_chr(disease_ids, function(disease_id) {
    if (num_failuers < get_cavalier_opt('hpo_api_max_failuers')) {
      result <- tryCatch(
        hpo_api_get(disease_id, prefix = 'network/annotation'),
        error = function(e) { 'ERROR'}
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

#' @export
build_disease_name_cache <- function(ORPHA = FALSE) {
  
  fun <- function() {
    disease_name_tbl <- 
      get_gene_disease_map(source = 'OMIM') %>% 
      select(disease_id) %>% 
      distinct()
    
    if (ORPHA) {
      disease_name_tbl <- 
        bind_rows(
          disease_name_tbl,
          get_gene_disease_map(source = 'ORPHA') %>% 
            select(disease_id) %>% 
            distinct()
        )
    }
    
    disease_name_tbl %>% 
      mutate(disease_name = hpo_api_get_disease_names(disease_id)) %>% 
      na.omit()
  }
  
  cache(
    fun = fun,
    name = 'disease_names',
    disk = TRUE,
    overwrite = TRUE,
    subdir = 'HPO'
  )
}

get_disease_names <- function(disease_ids) {
  
  tibble(disease_id = disease_ids) %>% 
    left_join(
      cache(
        fun = function() tibble(disease_id = character(), disease_name = character()), 
        name = 'disease_names',
        disk = TRUE, 
        store = FALSE
      ),
      by = 'disease_id'
    ) %>% 
    mutate(disease_name = replace(disease_name, 1:10, NA_character_)) %>% 
    (function(x) {
      missing <- which(is.na(x$disease_name)) 
      if (length(missing)) {
        x$disease_name[missing] <- hpo_api_get_disease_names(x$disease_id[missing])
      }
      x
    }) %>% 
    pull(disease_name)
}




